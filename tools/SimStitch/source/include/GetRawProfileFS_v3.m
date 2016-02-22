function [FS,FSParams,Instrument,NULL_REGION] = GetRawProfileFS_v3(fileList,options,file,DISPLAY_HDRS,X_ZFILLS)
%returns data and frequency points (increasing values) from .raw profile file

%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHANGE HISTORY
%version 1.1: removed criteria for all A parameters to be the same: now just issues a warning
%version 1.2: updated for new fileList structure
% 25/Oct/07: Added message to press any key on warning
%version 1.3
% 13/Mar/08:    updated for new f2mz and mz2f functions, now delete Raw COM
%               server
%
% Current Version 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%global constants

if isempty(X_ZFILLS), error('X_ZFILLS not defined'); end

%intialise
FS = [];
FSParams = [];

if ~isfield(options,'getSpec'), options.getSpec = 1; end
if ~isfield(options,'sumScans'), error('sumScans not specified'); end

%get Raw and Detector handles
%%%% hRaw = GetRawHandle(fileList.rawFull{file});
%%%% hDetector = GetDetectorHandle(hRaw);

if isunix 
	['wine C:\\python27\\python.exe ', which('MSFileReaderPy.py'), ' -i "' [fileList.RootDirectory, fileList.Samples(1,file).dataFile], '" -o "', fullfile(fileList.RootDirectory, 'temp.mat"')]
    system(['wine C:\\python27\\python.exe ', which('MSFileReaderPy.py'), ' -i "' [fileList.RootDirectory, fileList.Samples(1,file).dataFile], '" -o "', fullfile(fileList.RootDirectory, 'temp.mat"')]);
else
    system([which('MSFileReaderPy.exe'), ' -i "' [fileList.RootDirectory, fileList.Samples(1,file).dataFile], '" -o "', fullfile(fileList.RootDirectory, 'temp.mat"')]);
end

matPy = load(fullfile(fileList.RootDirectory, 'temp.mat'));

% Get instrument type
Instrument = matPy.temp.InstrumentName;
if strcmp(Instrument,'Orbitrap Elite')
    fileList.Instrument = 'orbitrap';
elseif strcmp(Instrument,'Q Exactive Orbitrap');
    fileList.Instrument = 'qexactive';
else
    fileList.Instrument = 'ltqft';
end
Instrument = fileList.Instrument; 

%get the list of scans for this spectrum
scanList = matPy.temp.ScanList;

%get handle to spectrum and header for each scan
%%%% hFilters = hDetector.get('Filters');
hFilters = matPy.temp.Filters;

%%%% filtersCount = hFilters.count;
filtersCount = matPy.temp.FiltersCount;

%if there is only one filter: easy
if filtersCount == 1
    filter = 1;
    %%%% hFilter = hFilters.Item(filter);
    hFilter = hFilters(1,:);
    %%%% if DISPLAY_HDRS, disp(['Filter: ', hFilter.Text]); end
    if DISPLAY_HDRS, disp(['Filter: ', hFilter]); end
    filterOfScanIndex = repmat(filter,1,length(scanList));
%otherwise need to find the filter index for each scan
else
    for i=1:length(scanList)       
        %%%% hFilter = hFilters.ScanNumber(scanList(i));
        hFilter = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(scanList(i)) ,'.Filter'));
        targetFilter = hFilter;
        %%%% targetFilter = hFilter.Text;
        %%%% if DISPLAY_HDRS, disp(['Filter: ', hFilter.Text,', scan: ',num2str(i)]); end
        if DISPLAY_HDRS, disp(['Filter: ', hFilter,', scan: ',num2str(i)]); end
        
        found = 0;
        for filterIndex = 1:filtersCount
            %%%% hFilter = hFilters.Item(filterIndex);
            hFilter = hFilters(filterIndex,:);
            [~,bracketpos] = find(hFilter == ']');
            if bracketpos < length(hFilter)
                hFilter(bracketpos+1:end) = [];
            end
            %%%% if isequal(hFilter.Text,targetFilter), found=1; break; end
            if isequal(char(hFilter),targetFilter), found=1; break; end
        end
        if ~found, error('Cannot find filter'); end
        filterOfScanIndex(i) = filterIndex;   
    end
end

% Determine m/z bounds and range of each filter
for i = 1:filtersCount
    [~, pos1] = find(hFilters(i,:) == '['); [~, pos2] = find(hFilters(i,:) == '-');
    mzlow(i,:) = str2num(hFilters(i,pos1+1:pos2(end)-1));
    [~, pos3] = find(hFilters(i,:) == ']');
    mzhigh(i,:) = str2num(hFilters(i,pos2(end)+1:pos3-1));
    mzrange(i,:) = mzhigh(i,:) - mzlow(i,:);
end

switch(fileList.Instrument)
    case{'ltqft'}
         filtersUsed = unique(filterOfScanIndex);
    case{'qexactive','orbitrap'}
        filtersUsed = unique(filterOfScanIndex); % Unique filters
        [~,filterpos] = find(filterOfScanIndex == 2);
        filterpos = filterpos(1)-1; % Detect last window spray stabilization time
        RepFilter = diff(filterOfScanIndex(filterpos(end)+1:end));% Detect repeating windows
        [~,RepFilter] = find(RepFilter < 1); 
        if ~isempty(RepFilter) %Windows are repeating
            RepFilter = RepFilter(1);
            filterOfScanIndex([filterpos,filterpos(end)+RepFilter+1:length(filterOfScanIndex)]) = 1; % Set all filters not to be used to index of first filter
        end
        rangedif = abs(mzrange(end) - mzrange); % Determine absolute difference mz range of last filter to all other filters
        if rangedif(1) < min(rangedif(2:end-1));
            filtersUsed(end) = [];
            mzlow(end) = []; mzhigh(end) = []; mzrange(end) = [];
        end
        clear rangedif;
        rangedif = abs(mzrange(2) - mzrange); % Determine absolute difference mz range of second filter to all other filters
        if rangedif(1) < min(rangedif(3:end)) % Check if second header is also a full scan
                filtersUsed([1,2]) = []; % Remove index first and second filter; scans with filterindex = 1 or 2 are not used anymore.
                mzlow([1,2],:) = []; mzhigh([1,2],:) = [];
        else
            filtersUsed(1) = []; % Remove index first filter; scans with filterindex = 1 are not used anymore.
            mzlow(1,:) = []; mzhigh(1,:) = []; 
        end
        %filterOfScanIndex([1:72,97:end]) = 1e12; % Temporary code to
        %process data Elliot (04/01/2016)
end

filtersCount = length(filtersUsed);

% Determine window overlap for stitching
NULL_REGION.MZMIN = 0; NULL_REGION.MZMAX = 0; % Dont cut off beginning for first window and end of last window
count = 1;
for i = 1:filtersCount
    if i ~= 1
        overlap(count) = (mzhigh(i-1) - mzlow(i))/2;
        bound(count) = (mzhigh(i-1) + mzlow(i-1))/2;
        count = count +1;
    end
end
mzbound = [mzlow(1,:)];
for i = 1:length(overlap)
    if i == 1
        mzstart = overlap(1);
        mzend = overlap(1);
    elseif overlap(i) ~= overlap(i-1)
        mzstart = [mzstart overlap(i)];
        mzend = [mzend overlap(i)];
        mzbound = [mzbound bound(i)];
    end
end
mzbound = [mzbound mzhigh(end,:)];
NULL_REGION.BOUND = mzbound;
NULL_REGION.START = mzstart;
NULL_REGION.END = mzend;

% STITCH.
spectrum.centroidsMZ = [];
spectrum.centroidsData = [];

%get the spectra
if DISPLAY_HDRS, disp(['File: ',fileList.Samples(1,i).ID,'...']); end
tic; lastT=0;

%loop through each filter (these will form the spectra)
spectraScans = {};
spectraFilter = {};

%options.sumScans = 1; (for QE same results should be obtained with
%options.sumScans = 0;)
if options.sumScans
    
    for k = 1:filtersCount

        filter = filtersUsed(k);
        scansThisFilter = [];
        d = [];

        %first get a list of scans from this filter
        idx = filterOfScanIndex==filtersUsed(k);
        filterScans = scanList(idx);
        scanIT = [];
        scanTIC = [];

        %get the handle to the spectra for this filter
        %%%% hSpectra = hDetector.get('Spectra', filter);

        for i=1:length(filterScans)
            
            %get the spectrum for this scan
            %first check for varying conversion parameters by extracting the trailer extra info
            switch(fileList.Instrument)
                case{'ltqft'}
                    newParams.A = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.A'))*1e3;% Multiply by 1e3 since equation is made for frequency in kHz while we use frequency in Hz
                    if i==1, fParams.A = newParams.A;
                    elseif newParams.A~=fParams.A
                        warning(['Parameter A is changing in ',fileList.Samples(1,file).ID,' - included, but check scan ',num2str(filterScans(i))]);
                    end
                    
                    newParams.B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6; 
                    if i==1, fParams.B = newParams.B;
                    elseif newParams.B~=fParams.B
                        warning(['Parameter B is changing in ',fileList.Samples(1,file).ID,' - included, but check scan ',num2str(filterScans(i))]);
                    end
                case{'qexactive','orbitrap'}
                    newParams.B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
                    if i==1, fParams.B = newParams.B;
                    elseif newParams.B~=fParams.B
                        warning(['Parameter B is changing in ',fileList.Samples(1,file).ID,' - included, but check scan ',num2str(filterScans(i))]);
                    end
                    
                    newParams.C = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.C'))*1e12; 
                    if i==1, fParams.C = newParams.C;
                    elseif newParams.C~=fParams.C
                        warning(['Parameter C is changing in ',fileList.Samples(1,file).ID,' - included, but check scan ',num2str(filterScans(i))]);
                    end
            end
            scanTIC = [scanTIC eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.TIC'))];
            
            scanIT = [scanIT eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'))];  
            
            switch(fileList.Instrument)
                case{'ltqft'}
                    if eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT')) >= 750 % Change for QExactive?.
                        warning(['Underfill in ',fileList.Samples(1,file).ID,' - NOT!!! ignoring scan ',num2str(filterScans(i))]);
                    end
            end
            
            if options.getSpec
                %get spectrum
                %can store the data values as singles to save space - can't do the same with m/z values, need the precision

                scanmz = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scanmz'));
                scand = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scand'));
                
                % Interpolation of noise packets values to m/z range of
                % interest
                if length(eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise'))) > 1% Interpolation when at least two noise packets are present
                    scanNoise = interp1(eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Mass')),eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise')),scanmz,'linear');
                    scanBase = interp1(eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Mass')),eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Base')),scanmz,'linear');
                    % Set value outside of bounds noise packets.
                    nan = ~isnan(scanNoise);
                    first_idx = find(nan == 1,1,'first');
                    last_idx = find(nan == 1,1,'last');
                    scanNoise(1:first_idx-1) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise(1)'));
                    scanNoise(last_idx+1:end) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise(end)'));
                    scanBase(1:first_idx-1) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Base(1)'));
                    scanBase(last_idx+1:end) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Base(end)'));
                elseif length(eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise'))) == 1
                    scanNoise = (eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise')))*ones(1,length(scanmz));
                    scanBase = (eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Base')))*ones(1,length(scanmz));
                else 
                    scanNoise = quantile(scand(scand~=0),0.01)*ones(1,length(scanmz)); % Estimate Noise value as 1% quantile of nonzero values in scand (the measured peak intensities)
                    scanBase = zeros(1,length(scanmz)); % scanBase estimated as zero
                end
                scanSNR = (scand-scanBase)./scanNoise;
                
                
                switch(fileList.Instrument)
                    case{'ltqft'}
                        %if the current A or B parameters have changed (are not the selected global parameters), normalise the m/z values to the global A and B parameter
                        if newParams.B~=fParams.B || newParams.A~=fParams.A 
                            scanmz = mz2f(scanmz,newParams,fileList.Instrument); 
                            scanmz = f2mz(scanmz,fParams,fileList.Instrument); 
                        end
                    case{'qexactive','orbitrap'}
                        if newParams.B~=fParams.B || newParams.C~=fParams.C 
                            scanmz = mz2f(scanmz,newParams,fileList.Instrument); 
                            scanmz = f2mz(scanmz,fParams,fileList.Instrument); 
                        end
                end
                             
                %store scan data
                if isempty(d)
                    d = scand;
                    mz = scanmz;
                    
                    % TEMP
                    noise = scanNoise;
                    base = scanBase;
                    SNR = scanSNR;
                else    %append data   %TO DO: WHAT TO DO WHEN I WANT TO COMBINE SCANS?
                    newd = [];
                    
                    % TEMP
                    newNoise = [];
                    newBase = [];
                    newSNR = [];
                    for j=1:size(d,1)
                        %in the case that there are multiple scans, we need to keep the
                        %number of data point entries the same, such that all d's can
                        %use the same 'mz'
                        [newd(j,:),scand_temp,newmz] = FillOut(d(j,:),mz,scand,scanmz); 
                        
                        % TEMP
                        [newNoise(j,:),scanNoise_temp,~] = FillOut(noise(j,:),mz,scanNoise,scanmz); 
                        [newBase(j,:),scanBase_temp,~] = FillOut(base(j,:),mz,scanBase,scanmz); 
                        [newSNR(j,:),scanSNR_temp,~] = FillOut(SNR(j,:),mz,scanSNR,scanmz); 

                    end
                    d = sum([newd;scand_temp],1);   %sum here to increase efficiency
                    mz = newmz;
                    
                    % TEMP
                    noise = sum([newNoise;scanNoise_temp],1);
                    base = sum([newBase;scanBase_temp],1);
                    SNR = sum([newSNR;scanSNR_temp],1);
                end
            end
            
            
            scansThisFilter = [scansThisFilter filterScans(i)];

        end

        spectraScans{k} = scansThisFilter;
        spectraFilter{k} = filter;

        %average the parameters if summing spectra over each filter: save each value otherwise
        IT(k) = mean(scanIT);
        TIC(k) = mean(scanTIC);
        % TEMP take average of noise, baseline and SNR! (replace with
        % normalization step + rescale to average of maximum peak?)
        
        %%%% mzStart(k) = hSpectrumHeader.LowMass;
        mzStart(k) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.LowMass'));
		
        %%%% mzEnd(k) = hSpectrumHeader.HighMass;
        mzEnd(k) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.HighMass'));
        
        params(k) = fParams;

        if options.getSpec
            %average spectra abundance
            data = d ./ length(scansThisFilter);
            %convert to frequency (increasing values); decalibrate
            spectrum(k).f = mz2f(mz,params(k),fileList.Instrument); %Update for other machines
            spectrum(k).f = fliplr(spectrum(k).f);
            spectrum(k).data = fliplr(data);
            spectrum(k).noise = fliplr(noise)./ length(scansThisFilter);
            spectrum(k).base = fliplr(base)./ length(scansThisFilter);
            spectrum(k).snr = fliplr(SNR)./ length(scansThisFilter);
        end

        %progress counter
        progress = k/filtersCount*100;
        t=toc;
        if t<5, lastT = floor(t/10);
        elseif t>5 && floor(t/10)>lastT
            disp([num2str(progress),'%']);lastT=floor(t/10);drawnow;
        end
        
    end
    numberOfSpectra = filtersCount;

else    %not summing scans TO DO OOK NOISE ETC AANPASSEN!!!

    scanCount = 0;
    for k = 1:filtersCount

        filter = filtersUsed(k);

        %first get a list of scans from this filter
        idx = filterOfScanIndex==filtersUsed(k);
        filterScans = scanList(idx);

        %get the handle to the spectra for this filter

        %loop through each scan in this filter
        for i=1:length(filterScans)

            %get the spectrum for this scan
            scanCount = scanCount + 1;
            nScans(scanCount) =1;
            spectraScans{scanCount} = filterScans(i);
            spectraFilter{scanCount} = filter;

            TIC(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.TIC'));
            IT(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'));
            mzStart(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.LowMass'));
            mzEnd(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.HighMass'));
            switch(fileList.Instrument)
                case{'ltqft'}
                    params(scanCount).A = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.A'))*1e3; % Multiply by 1e3 since equation is made for frequency in kHz while we use frequency in Hz
                    params(scanCount).B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
                case{'qexactive','orbitrap'}
                    params(scanCount).B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
                    params(scanCount).C = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.C'))*1e12;
            end
            if options.getSpec
                %get spectrum
                
                %can store the data values as singles to save space - can't do the same with m/z values, need the precision
                scanmz = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scanmz'));
                scand = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scand'));
                
                scanNoise = interp1(eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Mass')),eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Noise')),scanmz,'linear');
                scanBase = interp1(eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Mass')),eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.NoisePackets.Base')),scanmz,'linear');
                scanSNR = (scand-scanBase)./scanNoise;
                
                %progress counter
                progress = i/length(filterScans)*100;
                t=toc;
                if t<5, lastT = floor(t/10);
                elseif t>5 && floor(t/10)>lastT
                    disp([num2str(progress),'%']);lastT=floor(t/10);drawnow;
                end

                %convert to frequency (increasing values); decalibrate and store
                spectrum(scanCount).f = mz2f(scanmz,params(scanCount),fileList.Instrument); %Update for other machines
                spectrum(scanCount).f = fliplr(spectrum(scanCount).f);
                spectrum(scanCount).data = fliplr(scand);
                spectrum(scanCount).snr = scanSNR;
                spectrum(scanCount).noise = scanNoise;
                spectrum(scanCount).base = scanBase;
            end
        end
    end
    numberOfSpectra = scanCount;
end
    
%now pad out (uncompress) the data and frequency values
%theoretical delta f is BW/n, but we don't know BW/n
if options.getSpec
    for i=1:numberOfSpectra

        %first estimate the best frequency step to use (there are slight variations
        %in the actual frequency step probably due to representation and rounding errors.
        deltaf = median(spectrum(i).f(2:end)-spectrum(i).f(1:end-1));
        T(i) = 1/deltaf/2^X_ZFILLS; %This value doesnt really make sense for the QE since we do not not how many zero filles were carried out
        %if round(T(i)*1000) ~= 768, warning(['T not 0.768 in raw file, is ',num2str(T(i))]); end; %keyboard; end;

        %each data point will have a new index relative to the first frequency
        newIdx = (spectrum(i).f-spectrum(i).f(1))./deltaf;  %index relative to first freq. val
        % WARNING NOT RELEVANT FOR PEAK PICKED D
        %if find(abs(round(newIdx)-newIdx)>0.3), warning(['Possible error finding data indices in GetRawProfileFS.m! Press any key to continue...']); end %pause; end
        newIdx = round(newIdx) + 1; 
        
        %the new frequency values are the indices multiplied by the frequency step, plus the starting frequency
        spectrum(i).newf = ((0:max(newIdx)-1).*deltaf) + spectrum(i).f(1);
        
        % 19-08-2015: JE: Use interpolation to find new intensities etc. (The QE seems to use multiple frequency steps around. More accurate results when interpolation is used compared to old code)
        spectrum(i).newData = interp1(spectrum(i).f,spectrum(i).data,spectrum(i).newf);
        spectrum(i).newnoise = interp1(spectrum(i).f,spectrum(i).noise,spectrum(i).newf);
        spectrum(i).newbase = interp1(spectrum(i).f,spectrum(i).base,spectrum(i).newf);
        spectrum(i).newSNR = (spectrum(i).newData - spectrum(i).newbase) ./ spectrum(i).newnoise;

%         %spread out the data into the new data matrix, at the same time padding the gaps with zero
%         spectrum(i).newData(newIdx) = spectrum(i).data;
%         spectrum(i).newnoise(newIdx) = spectrum(i).noise;
%         spectrum(i).newbase(newIdx) = spectrum(i).base;
%         spectrum(i).newSNR(newIdx) = spectrum(i).snr;
% 
%         %the new frequency values are the indices multiplied by the frequency step, plus the starting frequency
%         spectrum(i).newf = ((0:length(spectrum(i).newData)-1).*deltaf) + spectrum(i).f(1);
        
    end
end

%return values of interest
for i=1:numberOfSpectra
    if options.getSpec
        FS(i).data = spectrum(i).newData;
        FS(i).f = spectrum(i).newf;
        FSParams(i).T = T(i);
        FS(i).Noise = spectrum(i).newnoise;
        FS(i).Base = spectrum(i).newbase;
        FS(i).SNR = spectrum(i).newSNR;
    end
    FSParams(i).mzStart = mzStart(i);
    FSParams(i).mzEnd = mzEnd(i);
    FSParams(i).TIC = TIC(i);
    switch(fileList.Instrument)
        case{'ltqft'}
            FSParams(i).A = params(i).A;
            FSParams(i).B = params(i).B;
            FSParams(i).C = 0;
        case{'qexactive','orbitrap'}
            FSParams(i).A = 0;
            FSParams(i).B = params(i).B;
            FSParams(i).C = params(i).C;
    end
    FSParams(i).zFills = X_ZFILLS;
    FSParams(i).IT = IT(i);
    FSParams(i).filterIndex = spectraFilter{i};
    FSParams(i).scans = spectraScans{i};
end
end

