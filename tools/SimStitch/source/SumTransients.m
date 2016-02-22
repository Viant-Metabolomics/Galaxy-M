function SumTransients(fileList, numscans, mintic, ignorecal, maxmz, html_outfile, avgd_trans_dir, outfile_messages)
%
%
% SumTransients_BB(fileList, numscans, mintic, ignorecal,maxmz, outfile)
%
% The 1st step in SimStitch processing of DIMS FTICR data.
%
% inputs:
%	fileList:   	the full path of an xml file containing information regarding location of the RAW files etc
%	numscans:   	an integer describing the number of scans to be summed. Set to 0 for 'all'
%	mintic:     	a percentage value describing a minimum threshold for a scan to be considered usable. Mintic is a percentage of the maximum TIC and is commonly around 40 i.e. the TIC of any scan must be at least 40% of the max TIC for all scans 
%	ignorecal:  	a boolean (1 or 0) that is used to declare whether to ignore scans with varying A and B parameters 
%	maxmz:      	this value can be used to set the maximum m/z to be included. Set to 0 for automatic estimation.
%	html_outfile:  	full path for a file to be used by Galaxy for presenting links to the output files
%	avgd_trans_dir:	full path for the directory (defined by Galaxy) that should hold the output files referenced by html_file
%	outfile_messages: full path to a text file to store messages related to used/unused transients etc (available in Galaxy history)
%	outfile_params:	full path to a text file to store the params, making them available in Galaxy history.  
%
%
% outputs:
%	this function does not return any Matlab variables because it is designed to create various output files as expected by Galaxy. 
%
%
%
% This file has been adapted for Galaxy from many previous incarnations by R.L.Davidson
% For more information contact Mark Viant m.viant@bham.ac.uk 
%
%
% Possible TO ADD:
%	1) output the unused transients to an HTML file and directory
%	2) clean up the code! remove legacy code.
%	3) reduce fileList to only the common variables that are required for this script. Now that input/output files are stored by Galaxy there is no need for all the directory details in fileList.
%
%	 

if isa(numscans, 'char')
	numscans = str2num(numscans); %number of scans to read in (set to 0 for all scans)
end

if isa(mintic, 'char')
	mintic = str2num(mintic); %the minimum allowable TIC for a scan to be included (as % of maximum TIC), typically ~40
end

if isa(ignorecal, 'char')
	ignorecal = str2num(ignorecal);
end

if isa(maxmz, 'char')
	maxmz = str2num(maxmz); %manually set the maximum m/z (set to 0 for automatic max m/z taken from the .raw file)
end

%extract the fileList information from the xml file
fileList = ImportFileListXML(fileList);

%set some basic parameters required later. These are never changed in regular practice so were not added to Galaxy input form.
DISPLAY_HDRS = 0;  %display file header info
X_ZFILLS = 1;      %zero-fills assumed to have been applied by xCalibur for .raw file spectra

%counter for output file naming
segCount = 0;
batchMsg = {};

fileList.ATDir = [avgd_trans_dir,filesep] ; %this is for Galaxy output

try mkdir(fileList.ATDir);
catch
    error(['Cannot make averaged transient directory ',fileList.ATDir]); 
end

dFold = [fileList.ATDir,'unusedtransients'];
try mkdir(dFold);
catch
    error(['Cannot make unused transient directory ',dFold]); 
end

filename_html = {}; %this is for storing Galaxy output

%first get the parameter information ONLY for each scan from the Data file
for i=1:fileList.nDataFiles

    segCount = 0;
    message = {};

    % If else apparaat
    % read in raw file
    fprintf('\nFile: %s\n',fileList.Samples(1,i).ID);
    options.sumScans = 0;
    options.getSpec = 0;
    
    % Switch between .raw and .ser
    if strcmp(fileList.Samples(1,1).dataFile(end-2:end),'raw') || strcmp(fileList.Samples(1,1).dataFile(end-2:end),'RAW') % ltqft
            [~,specParamsR,Instrument,NULL_REGION] = GetRawProfileFS_v3(fileList,options,i,DISPLAY_HDRS,X_ZFILLS);
            fileList.Instrument = Instrument;
    else % Solarix
            % To do: Zoek method folder naam.
            specParamsR = Read_method_file([fileList.RootDirectory, fileList.Samples(1,i).dataFile,'']);
            fileList.Instrument = 'solarix';
    end

    %now loop through the different m/z ranges
    listFilterIndices = unique([specParamsR.filterIndex]);
    for f=1:length(listFilterIndices)
        
        %parameters for this SIM window
        specParams = specParamsR([specParamsR.filterIndex]==listFilterIndices(f));

        %check SIM window is required
        if maxmz && (specParams(1).mzEnd > maxmz)
            %SIM window not required
            message{end+1} = [fileList.Samples(1,i).ID,'_seg',num2str(listFilterIndices(f),'%.2d'),': beyond required range'];
        else
            segCount = segCount + 1;

            %now decide which of these scans we wish to read in
            numSpectra = length(specParams);
            scanList = 1:numSpectra;    %index
            maxTIC = max([specParams.TIC]);
            fprintf('Segment %d ignoring scans:',f);
            for j=1:numSpectra
                if specParams(j).TIC/maxTIC*100 < mintic
                    ulf=specParams(j).TIC/maxTIC*100
                    fprintf('%d\t',specParams(scanList(j)).scans);
                    message{end+1} = [fileList.Samples(1,i).ID,'_seg',num2str(listFilterIndices(f),'%.2d'),': IGNORING SCAN ',num2str(specParams(scanList(j)).scans),' (TIC is <40% max seg TIC: RAW TIC = ',num2str(specParams(j).TIC),', MAX TIC = ',num2str(maxTIC),')'];
                    % copy all poor scans to a folder
                    dFold = [fileList.ATDir,'unusedtrans'];
                    if ~isdir(dFold)
                        [status,~,~] = mkdir(dFold);
                        if ~status, error('Cannot create folder for transients'); end
                    end
                    sFile = [[fileList.Samples(1,i).ID,'_'],num2str(specParams(scanList(j)).scans),'.dat'];
                    fprintf(' - copying %s to %s\n',sFile,dFold);
                    [status,~,~] = copyfile([fileList.RootDirectory,sFile],dFold);
                    if ~status, error('Cannot copy transient'); end
                    scanList(j) = 0;
                end
            end
            if isempty(find(scanList==0, 1)), fprintf(' (none)');
            else scanList = scanList(scanList~=0);    %relative
            end
            fprintf('\n');
            scanList = [specParams(scanList).scans];   %absolute
            message{end+1} = [fileList.Samples(1,i).ID,'_seg',num2str(listFilterIndices(f),'%.2d'),': ',num2str(length(scanList)),' scans: ',num2str(scanList)];
            
            %now read in the transient files for this segment
            scanListAll = scanList;
            if numscans && (length(scanList) > numscans)
                scanList = scanList(1:numscans);
                warning(['First ',num2str(numscans),' (compliant) scans only']);
            end
            if strcmp(fileList.Samples(1,1).dataFile(end-2:end),'raw') || strcmp(fileList.Samples(1,1).dataFile(end-2:end),'RAW')% ltqft
                [transient,specParamsT,segMessage] = ReadInDat_3_0(fileList.RootDirectory,[fileList.Samples(1,i).ID,'_'],scanList,'all',ignorecal,DISPLAY_HDRS);
            else
                % Read in solarix .ser files
            end
            specParamsT.NULL_REGION  = NULL_REGION; % Store window overlap info for use in stitch.m
            if isempty(transient), error('empty transient'); end
            message = [message segMessage];
            %scans that should have been read in but failed
            scansNotRead = setdiff(specParamsT.scans,specParamsT.scansRead);
            %copy out
            for s=1:length(scansNotRead)
                sFile = [fileList.RootDirectory,[fileList.Samples(1,i).ID,'_'],num2str(scansNotRead(s)),'.dat'];
                fprintf(' - copying %s to %s\n',sFile,dFold);
                [status,~,~] = copyfile(sFile,dFold);
                if ~status, error('Cannot copy transient'); end
            end            

            % copy unrequired scans out
            for s=1:length(scanListAll)
                if isempty(find(scanList==scanListAll(s), 1))
                    sFile = [fileList.RootDir,[fileList.Samples(1,i).ID,'_'],num2str(specParams(scanListAll(s)).scans),'.dat'];
                    fprintf(' - copying %s to %s\n',sFile,dFold);
                    [status,~,~] = copyfile(sFile,dFold);
                    if ~status, error('Cannot copy transient'); end
                end
            end

            %merge some parameters from the raw file
            specParamsT.scans = [specParams.scans];
            specParamsT.TICraw = [specParams.TIC];
            specParamsT.IT = [specParams.IT];
            specParamsT.filterIndex = specParams(1).filterIndex;
            specParamsT.mzStart = specParams(1).mzStart;
            specParamsT.mzEnd = specParams(1).mzEnd;

            %save transient
            fileName = [fileList.ATDir,fileList.Samples(1,i).ID,'_seg',num2str(segCount,'%.2d')];
            fprintf('\tSaving averaged transients to %s...',fileName);
            save(fileName,'transient','specParamsT');
            fprintf('done\n');

	    %prepare to save to html_out file for Galaxy
	    filename_html{end+1} = [fileList.Samples(1,i).ID,'_seg',num2str(segCount,'%.2d'), '.mat'];


        end %check SIM window required
    end
    %save messages (UPDATED FOR GALAXY)
    fprintf('Saving message file: %s\r\n',outfile_messages);
    fid = fopen(outfile_messages,'a');
    if ~fid, error('Cannot create message file'); end
    for m=1:length(message)
        fc=fprintf(fid,'%s\r\n',message{m});
        if ~fc, error('Cannot write to message file'); end
    end
    fclose(fid);
end

%write html for Galaxy
fid = fopen(html_outfile, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Sum Transient Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>Sum Transient Data - Output</b></h2></div>');
fprintf(fid,'<hr></hr>')
for i=1:length(filename_html)
  	html_code = ['<div style="padding: 3px"><b><a href="',filename_html{i},'">',filename_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);


