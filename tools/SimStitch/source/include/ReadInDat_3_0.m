function [transient,specParams,message] = ReadInDat_3_0(fDir,fStem,scanList,n,ignorecal,DISPLAY_HDRS)
%read in n data points from a list of .dat transient files in scanList
% 
%given:
% fStem file stem
% n, number of data points to read in (n = 'all' -> read all data available)
% scanlist: if scanList is empty, read in all available transients
% if ignorecal=1, transients with different
% calibration parameters from the first scan in the segment are included,
% and the output calibration parameters are the average of the parameters
% from the transients, otherwise they are not included
%
%outputs:
% transient = row vector length n
% specParams = BW VT A B C TIC T scans and scansread: for each spectrum (scans is list
    %of all 'potential' scan numbers, scansread is list of scan numbers
    %actually summed into transient)
% message = text message giving information regarding scans (not) read

%UPDATES:
%29/May/07 v0.1: Added list of *potential* scans (ie including any failed)
%6/Aug/07   v0.2: Removed options and changed storage of transient to
%single (32 bit) to save storage space. Matches 32 bit storage in dat
%files, so no precision loss
%5/Nov/07   v0.3: Added inherited parent parameters (for ignoring changing
%parameters) input, and scansNotRead and message outputs. Moved scansread
%to specParams.  Added check for non-zero C parameter.
%23/Nov/07  v0.4: Increased maximum transient value to 5000.
%01/Apr/08  v0.5: Fixed minor bug if file not found or scanList empty.
%                 Removed QUIET
%21/Aug/08        Minor changes to error messages

%globals
DISPLAY_HDRS;

transient = [];
scansRead = [];      %list of valid scans read-in
scanList = unique(scanList);  %sort ascending
scans = []; %list of all potential scans
specParams = [];
A = [];
B = [];
C = [];
message = {};
i = 1;  % current index of scanList
done = 0;   % 1 when all (or n) transients read in

% first scan number
if isempty(scanList)
    currScan = 1;
else
    currScan = scanList(1);
end

while ~done

    %try to read in next file
    fileName = [fDir,fStem,num2str(currScan),'.dat'];
    fprintf('\tReading in file %s\n',fileName);
    fid = fopen(fileName,'r','n');
    if (fid == -1)
        %failed
        if isempty(scanList)
            disp(['No more files: ',fStem,num2str(scansRead(end)),'.dat last in series']);
            done = 1;
        else
            warning('File %s not found! Press any key to continue...',fileName); %pause;
            message{end+1} = ['File ',fileName,' not found, skipped'];
        end
    else
        includeThisScan = 1;
        
        %potential scan
        scans = [scans currScan];

        %check parameters
        %read in header info
        fseek(fid,0,'bof');
        Astr = [];
        for j=1:3
            An = 0;
            while An ~= 10
                An=fread(fid,1);
                Astr = [Astr An];
            end
        end
        if DISPLAY_HDRS, char(Astr); end
        %read parameters
        Bstr = cell(12,2);
        for j=1:12
            Bn = -1;
            Bchar = [];
            while Bn ~= double(':')
                Bn = fread(fid,1);
                Bchar = [Bchar Bn];
            end
            Bstr{j,1} = char(Bchar);
            Bn = -1;
            Bchar = [];
            while Bn ~= 10
                Bn = fread(fid,1,'*char');
                Bchar = [Bchar Bn];
            end
            Bstr{j,2} = Bchar(1:end-1);
        end
        if DISPLAY_HDRS
            for j = 1:size(Bstr,1)
                disp([Bstr{j,1},Bstr{j,2}]);
            end
        end
        %extract some parameters
        newBW = str2double(Bstr{5,2});     %bandwidth
        newVT = str2double(Bstr{11,2});    %trapping voltage
        maxN = str2double(Bstr{4,2});      %number of data points
        if (strcmp(n,'all') || (n < 0))
            newN = maxN;
        elseif (n > maxN)
            warning(['Maximum number of input data points exceeded (asked for ',num2str(n),')! Press any key to continue...']); %pause;
            message{end+1} = ['File ',fileName,' maximum number of data points exceeded, set to ',num2str(maxN)];
            newN = maxN;
        else
            newN = n;
        end
        newA = str2double(Bstr{8,2});      %conversion parameter A
        newB = str2double(Bstr{9,2});      %conversion parameter B
        newC = str2double(Bstr{10,2});     %conversion parameter C
        % not using C, so just check it's zero
        if newC, close(fid); error('C parameter non-zero in transient file'); end
        % check BW, VT and n haven't changed
        if (exist('BW') && (BW~=newBW)) || (exist('VT') && (VT~=newVT)) || (exist('N') && (N~=newN))
            close(fid); error('Essential transient parameters changing!');
        end
        % look for changing calibration parameters in transients
        if (~isempty(A) && (A(end)~=newA)) || (~isempty(B) && (B(end)~=newB)) || (~isempty(C) && (C(end)~=newC))
            warning(['Parameters changed in scan ',num2str(currScan),', possible underfill']);
            fprintf('A: %f to %f\n',A(end),newA);
            fprintf('B: %f to %f\n',B(end),newB);
            fprintf('C: %f to %f\n',C(end),newC);
            if ignorecal
                includeThisScan = 1;
                message{end+1} = ['File ',fileName,' A/B/C parameter change: included'];
            else
                includeThisScan = 0;
                fprintf('Ignoring scan. Press any key to continue...\n'); %pause;
                message{end+1} = ['File ',fileName,' A/B/C parameter change: NOT included'];
            end
        end

        % error check transient
        if includeThisScan
            %read in data
            if DISPLAY_HDRS, disp(['Reading in ',num2str(n), ' data points...']); end
            transientF = (single(fread(fid, newN, 'float32')))';
            %error check length
            if (length(transientF) ~= newN)
                warning([fileName,': mismatch in length! Ignoring file! Press any key to continue...']);%pause;
                message{end+1} = ['File ',fileName,' mismatch in length: NOT included'];
                includeThisScan = 0;
            %check for corrupt data
            elseif (max(abs(transientF)) > 5e3)
                warning([fileName,': transient values > 5000, maybe corrupt! Ignoring file! Press any key to continue...']);%pause;
                message{end+1} = ['File ',fileName,' transient values > 5000: NOT included'];
                min(transientF)
                max(transientF)
                includeThisScan = 0;
            end
        end

        % average transient data
        if includeThisScan
            if isempty(transient), transient = transientF;
            else transient = transient + transientF; %sum for averaging at end
            end
            scansRead = [scansRead currScan];
            scanCount = length(scansRead);
            %total ion current
            lengthn = length(transientF);
            TIC(scanCount) = norm(transientF)/sqrt(lengthn);
            %params that may vary
            A(scanCount) = newA;
            B(scanCount) = newB;
            C(scanCount) = newC;
            %supposedly constant parameters
            if ~exist('BW'), BW = newBW; end
            if ~exist('VT'), VT = newVT; end
            if ~exist('N'), N = newN; end
        end

        % done
        fclose(fid);
    end
    
    if ~done
        % next scan
        if isempty(scanList)
            currScan = currScan+1;
        else
            i = i+1;
            if i<=length(scanList)
                currScan = scanList(i);
            else
                done = 1;
            end
        end
    end
end

scanCount = length(scansRead);
fprintf('\tRead in %d scans\n',scanCount);

%average scans
transient = transient./scanCount;

%save out the parameters in structure
specParams.BW = BW;
specParams.VT = VT;
specParams.A = mean(A);
specParams.B = mean(B);
specParams.C = mean(C);
specParams.TIC = TIC;
ts = 1/(2*BW);
specParams.T = ts*N;        %NOTE T should include the lead-in time for the first sample!
                            %which is why T=ts*n and not T=ts*(n-1)
specParams.scans = scans;   %all potential scans
specParams.scansRead = scansRead;

if isempty(message), message = {'OK'}; end

return
