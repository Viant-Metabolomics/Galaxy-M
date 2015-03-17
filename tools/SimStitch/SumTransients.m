function SumTransients_BB(fileList, numscans, mintic, ignorecal, maxmz, html_outfile, avgd_trans_dir, outfile_messages)
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
%	html_file:  	full path for a file to be used by Galaxy for presenting links to the output files
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


%set input parameters to fit existing legacy code
SET = fileList; %this needs CHANGED now that an individual xml file is being passed

PARAMS.NSCANS = numscans; %number of scans to read in (set to 0 for all scans)
PARAMS.TICTHRESH_PC = mintic; %the minimum allowable TIC for a scan to be included (as % of maximum TIC), typically ~40
PARAMS.IGNORECHANGINGCALPARAMS_ON = ignorecal; %ignore scans which have varying A and B parameters? (Usually caused by underfill)
PARAMS.FIXEDMAX_MZ = maxmz; %manually set the maximum m/z (set to 0 for automatic max m/z taken from the .raw file)

%set some basic parameters required later. These are never changed in regular practice so were not added to Galaxy input form.
DISPLAY_HDRS = 0;  %display file header info
X_ZFILLS = 1;      %zero-fills assumed to have been applied by xCalibur for .raw file spectra


%extract the fileList information from the xml file
fileList = import_filelist_xml(SET);

%rawCount = length(fileList.rawFull);
%if fileList.rawContig
%    fprintf('Summing transients for %d contiguous RAW files...',rawCount);
%else
%    fprintf('Summing transients for %d separate RAW files...',rawCount);
%end

%counter for output file naming
segCount = 0;
batchMsg = {};

%make directory for summed transients and unused transients


fileList.avTransDir = [avgd_trans_dir,filesep] ; %this is for Galaxy output
dFold = fileList.avTransDir;
try mkdir(dFold);
catch
    error(['Cannot make averaged transient directory ',dFold]); 
end

dFold = [fileList.avTransDir,'unusedtrans_',fileList.setUID];
try mkdir(dFold);
catch
    error(['Cannot make unused transient directory ',dFold]); 
end

rawCount = length(fileList.rawFull);

filename_html = {}; %this is for storing Galaxy output

%first get the parameter information ONLY for each scan from the RAW file
for i=1:rawCount

    segCount = 0;
    message = {};

    % read in raw file
    fprintf('\nFile: %s\n',fileList.rawShort{i});
    options.sumScans = 0;
    options.getSpec = 0;

    [spec,specParamsR] = GetRawProfileFS_v2(fileList,options,i,DISPLAY_HDRS,X_ZFILLS);

%now loop through the different m/z ranges
    listFilterIndices = unique([specParamsR.filterIndex]);
    for f=1:length(listFilterIndices)
        
        %parameters for this SIM window
        specParams = specParamsR(find([specParamsR.filterIndex]==listFilterIndices(f)));

        %check SIM window is required
        if PARAMS.FIXEDMAX_MZ & (specParams(1).mzEnd > PARAMS.FIXEDMAX_MZ)
            %SIM window not required
            message{end+1} = [fileList.rawShort{i},'_seg',num2str(listFilterIndices(f),'%.2d'),': beyond required range'];
        else
            segCount = segCount + 1;

            %now decide which of these scans we wish to read in
            numSpectra = length(specParams);
            scanList = 1:numSpectra;    %index
            maxTIC = max([specParams.TIC]);
            fprintf('Segment %d ignoring scans:',f);
            for j=1:numSpectra
                if specParams(j).TIC/maxTIC*100 < PARAMS.TICTHRESH_PC
                    ulf=specParams(j).TIC/maxTIC*100
                    fprintf('%d\t',specParams(scanList(j)).scans);
                    message{end+1} = [fileList.rawShort{i},'_seg',num2str(listFilterIndices(f),'%.2d'),': IGNORING SCAN ',num2str(specParams(scanList(j)).scans),' (TIC is <40% max seg TIC: RAW TIC = ',num2str(specParams(j).TIC),', MAX TIC = ',num2str(maxTIC),')'];
                    % copy all poor scans to a folder
                    dFold = [fileList.avTransDir,'unusedtrans_',fileList.setUID];
                    if ~isdir(dFold)
                        [status,msg,msgid] = mkdir(dFold);
                        if ~status, error('Cannot create folder for transients'); end
                    end
                    sFile = [fileList.datStem{i},num2str(specParams(scanList(j)).scans),'.dat'];
                    fprintf(' - copying %s to %s\n',sFile,dFold);
                    [status,msg,msgid] = copyfile([fileList.datDir,sFile],dFold);
                    if ~status, error('Cannot copy transient'); end
                    scanList(j) = 0;
                end
            end
            if isempty(find(scanList==0)), fprintf(' (none)');
            else scanList = scanList(find(scanList));    %relative
            end
            fprintf('\n');
%             plot([specParams.scans],[specParams.TIC]);pause;
            scanList = [specParams(scanList).scans];   %absolute
            message{end+1} = [fileList.rawShort{i},'_seg',num2str(listFilterIndices(f),'%.2d'),': ',num2str(length(scanList)),' scans: ',num2str(scanList)];

            %now read in the transient files for this segment
            scanListAll = scanList;
            if PARAMS.NSCANS & (length(scanList) > PARAMS.NSCANS)
                scanList = scanList(1:PARAMS.NSCANS);
                warning(['First ',num2str(PARAMS.NSCANS),' (compliant) scans only']);
    %             batchMsg{end+1} = ['Reading ',num2str(length(scanList)),' scans in: ', fileList.datStem{i}];
            else
    %             batchMsg{end+1} = ['*Reading ',num2str(length(scanList)),' scans in: ', fileList.datStem{i}];
            end
            [transient,specParamsT,segMessage] = ReadInDat_3_0(fileList.datDir,fileList.datStem{i},scanList,'all',PARAMS, DISPLAY_HDRS);
            if isempty(transient), error('empty transient'); end
            message = [message segMessage];
            %scans that should have been read in but failed
            scansNotRead = setdiff(specParamsT.scans,specParamsT.scansRead);
            %copy out
            for s=1:length(scansNotRead)
                sFile = [fileList.datDir,fileList.datStem{i},num2str(scansNotRead(s)),'.dat'];
                fprintf(' - copying %s to %s\n',sFile,dFold);
                [status,msg,msgid] = copyfile(sFile,dFold);
                if ~status, error('Cannot copy transient'); end
            end            

            % copy unrequired scans out
            for s=1:length(scanListAll)
                if isempty(find(scanList==scanListAll(s)))
                    sFile = [fileList.datDir,fileList.datStem{i},num2str(specParams(scanListAll(s)).scans),'.dat'];
                    fprintf(' - copying %s to %s\n',sFile,dFold);
                    [status,msg,msgid] = copyfile(sFile,dFold);
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
            fileName = [fileList.avTransDir,fileList.rawShort{i},'_seg',num2str(segCount,'%.2d')];
            fprintf('\tSaving averaged transients to %s...',fileName);
            save(fileName,'transient','specParamsT');
            fprintf('done\n');

	    %prepare to save to html_out file for Galaxy
	    filename_html{end+1} = [fileList.rawShort{i},'_seg',num2str(segCount,'%.2d'), '.mat'];


        end %check SIM window required
    end
    %save messages (UPDATED FOR GALAXY)
    %fileName = [fileList.avTransDir,'avg_messages_',fileList.rawShort{i},'.txt'];
    fprintf('Saving message file: %s\r\n',outfile_messages);
    fid = fopen(outfile_messages,'w');
    if ~fid, error('Cannot create message file'); end
    for m=1:length(message)
        fc=fprintf(fid,'%s\r\n',message{m});
        if ~fc, error('Cannot write to message file'); end
    end
    fclose(fid);
end

%save parameters (probably not necessary when Galaxy is being used as the history and parameters can be saved from there)
%fileName = [fileList.avTransDir,'params.txt'];
%fprintf('Saving parameters to file: %s\n',outfile_params);
%fid = fopen(outfile_params,'w');
%if ~fid, error('Cannot create message file'); end
%fprintf(fid,'FILE VERSION:\t%s',mfilename);
%fprintf(fid,'\r\n\r\nPARAMETERS:\r\n');
%fn = fieldnames(PARAMS);
%for i=1:length(fn)
%    val = getfield(PARAMS,fn{i});
%    if isstr(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
%fclose(fid);
%disp('...finished');

% fprintf('Saving batch messages to %s\n',[fileList.avTransDir,'batchMsg']);
% save([fileList.avTransDir,'batchMsg'],'batchMsg');

%write html for Galaxy
fid = fopen(html_outfile, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Sum Transient Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
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
fprintf(fid,'<p>');
fprintf(fid,'Note: ---');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);
return



function fileList_struct = import_filelist_xml(file_location)
% 
% function fileList_struct = import_filelist_xml(file_location)
%
% function to allow SimStitch pipeline to read/import xml representations of individual fileList entries 
%
% specifically aimed at the Galaxy Project implementation of SimStitch etc.
% inputs:
%   file_location:  the location of a fileList.xml file as produced by
%                   FileListManagerGUI
%
% outputs:
%   fileList_struct:a matlab struct conforming to SimStitch's fileList
%                   structure
%
% R.L.Davidson 
% 29/05/2014



try
    tree = xmlread(file_location);
catch
    errordlg('Failed to read XML file.','Import error');
end

if ~tree.hasChildNodes()
    errordlg('XML tree empty!','Import error');
    return
end

% this struct declaration is not very concise, but helps me remember the structure (could be
% declared like fileListStruct.spec is e.g.
% struct('setIdentifier',{},'rootIdentifier',{},...)
fileListStruct = [];
fileListStruct.setIdentifier = '';
fileListStruct.rootDirectory= '';
fileListStruct.rawDirectory= [];
fileListStruct.spec= struct('ID',{},'rawFile',{},'sampleID',{});
fileListStruct.analysisSplit= [];
fileListStruct.datDirectory= [];
fileListStruct.avgdTransDirectory= '';
fileListStruct.overlappingSpecDirectory= '';
fileListStruct.stitchedSpecDirectory= '';
fileListStruct.peakListsDirectory= '';
fileListStruct.filteredPeaksDirectory= '';
fileListStruct.blankFlaggedDirectory= '';
fileListStruct.combinedPeaksDirectory= '';
fileListStruct.internalCalFile= [];
fileListStruct.numReps= [];
fileListStruct.count= [];
fileListStruct.notes= [];
fileListStruct.spectraCount= [];
fileListStruct.DSO_Directory= '';
fileListStruct.MultivarModel_Directory= '';
fileListStruct.MetID_Directory= '';
fileListStruct.specHeaders= {};
fileListStruct.specData= {};


%i mislabeled the struct above, will rename for ease now (lazy!)
fileList_struct = fileListStruct;

%start to interrogate the xml input
if strcmp('fileList',tree.getChildNodes.item(0).getNodeName)
    fileList = tree.getChildNodes.item(0);
else
    errordlg('XML Document root is not fileList - stopping')
    return
end

try
    
    if ~isempty(fileList.getElementsByTagName('setIdentifier').item(0).getChildNodes.item(0))
	fileList_struct.setIdentifier = char(fileList.getElementsByTagName('setIdentifier').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('rootDirectory').item(0).getChildNodes.item(0))    
	fileList_struct.rootDirectory = char(fileList.getElementsByTagName('rootDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('rawDirectory').item(0).getChildNodes.item(0))
        fileList_struct.rawDirectory = char(fileList.getElementsByTagName('rawDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    fileList_struct.analysisSplit = str2num(fileList.getElementsByTagName('analysisSplit').item(0).getChildNodes.item(0).getData);
    
    if ~isempty(fileList.getElementsByTagName('datDirectory').item(0).getChildNodes.item(0))
        fileList_struct.datDirectory = char(fileList.getElementsByTagName('datDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('avgdTransDirectory').item(0).getChildNodes.item(0))
        fileList_struct.avgdTransDirectory = char(fileList.getElementsByTagName('avgdTransDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('overlappingSpecDirectory').item(0).getChildNodes.item(0))
        fileList_struct.overlappingSpecDirectory = char(fileList.getElementsByTagName('overlappingSpecDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('stitchedSpecDirectory').item(0).getChildNodes.item(0))
        fileList_struct.stitchedSpecDirectory = char(fileList.getElementsByTagName('stitchedSpecDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('peakListsDirectory').item(0).getChildNodes.item(0))
        fileList_struct.peakListsDirectory = char(fileList.getElementsByTagName('peakListsDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('filteredPeaksDirectory').item(0).getChildNodes.item(0))
        fileList_struct.filteredPeaksDirectory = char(fileList.getElementsByTagName('filteredPeaksDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('blankFlaggedDirectory').item(0).getChildNodes.item(0))
        fileList_struct.blankFlaggedDirectory = char(fileList.getElementsByTagName('blankFlaggedDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('combinedPeaksDirectory').item(0).getChildNodes.item(0))
        fileList_struct.combinedPeaksDirectory = char(fileList.getElementsByTagName('combinedPeaksDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('internalCalFile').item(0).getChildNodes.item(0))
        fileList_struct.internalCalFile = char(fileList.getElementsByTagName('internalCalFile').item(0).getChildNodes.item(0).getData);
    end
    
    
    fileList_struct.numReps = str2num(fileList.getElementsByTagName('numReps').item(0).getChildNodes.item(0).getData);
    
    fileList_struct.count = str2num(fileList.getElementsByTagName('count').item(0).getChildNodes.item(0).getData);
    
    if ~isempty(fileList.getElementsByTagName('notes').item(0).getChildNodes.item(0))
        fileList_struct.notes = char(fileList.getElementsByTagName('notes').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('DSO_Directory').item(0).getChildNodes.item(0))
        fileList_struct.DSO_Directory = char(fileList.getElementsByTagName('DSO_Directory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('MetID_Directory').item(0).getChildNodes.item(0))
        fileList_struct.MetID_Directory = char(fileList.getElementsByTagName('MetID_Directory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('MultivarModel_Directory').item(0).getChildNodes.item(0))
        fileList_struct.MultivarModel_Directory = char(fileList.getElementsByTagName('MultivarModel_Directory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('specHeaders').item(0).getChildNodes.item(0))
        header_char = char(fileList.getElementsByTagName('specHeaders').item(0).getChildNodes.item(0).getData);
        specHeaders = {};
        while ~isempty(header_char)
            [specHeaders{end+1},header_char] = strtok(header_char,',');
        end
        fileList_struct.specHeaders = specHeaders;
    end
    
    if ~isempty(fileList.getElementsByTagName('specData').item(0).getChildNodes.item(0))
       data_char = char(fileList.getElementsByTagName('specData').item(0).getChildNodes.item(0).getData);
    
       data_char_cell = {};
       while ~isempty(data_char)
           [data_char_cell{end+1},data_char] = strtok(data_char,';');
       end
       
       num_cols = length(strfind(data_char_cell{1},','))+1;
       data_cell = cell(length(data_char_cell),num_cols);
       for i = 1:length(data_char_cell)
           data_char_row = data_char_cell{i};
           data_cell_row = {};
           while ~isempty(data_char_row)
               [data_cell_row{end+1},data_char_row] = strtok(data_char_row,',');
           end
           data_cell(i,:) = deal(data_cell_row);
       end
       
       dc = cellfun(@str2num,data_cell, 'UniformOutput', false);
       mask = cellfun(@isempty,dc);
       data_cell(~mask) = dc(~mask);
       
       fileList_struct.specData = data_cell;
    
    end
       
    
catch
    errordlg('XML file does not contain all fields - stopping!','Import error')
    return
end

try
    spec_node = fileList.getElementsByTagName('spec').item(0);
    instances = spec_node.getElementsByTagName('instance');
    
    for j = 1:instances.getLength;
        fileList_struct.spec(j).rawFile = char(instances.item(j-1).getElementsByTagName('rawFile').item(0).getChildNodes.item(0).getData);
        fileList_struct.spec(j).ID = char(instances.item(j-1).getElementsByTagName('ID').item(0).getChildNodes.item(0).getData);
        if ~isempty(instances.item(j-1).getElementsByTagName('sampleID').item(0).getChildNodes.item(0))
            fileList_struct.spec(j).sampleID = char(instances.item(j-1).getElementsByTagName('sampleID').item(0).getChildNodes.item(0).getData);
        end
    end
    
catch
    errordlg('XML fileList does not contain information about spectral raw files... quitting','Import error')
    return
end

%add some conversions to help this xml format fit with the expected
%fileList structure variable names from legacy code:

%SumTransients conversions
fileList_struct.avTransDir = fileList_struct.avgdTransDirectory;
fileList_struct.setUID = fileList_struct.setIdentifier;

fileList_struct.rawShort = {fileList_struct.spec.ID};

for i=1:length(fileList_struct.rawShort)
    fileList_struct.rawFull{i} = [fileList_struct.rootDirectory,fileList_struct.rawDirectory, fileList_struct.spec(i).rawFile];
end
fileList_struct.rawContig = 0; % set to 1 if raw files all combine to form a single spectrum (never used)

fileList_struct.datDir = [fileList_struct.rootDirectory, fileList_struct.datDirectory];
for i=1:length(fileList_struct.rawShort)
    fileList_struct.datStem{i} = [fileList_struct.rawShort{i},'_'];
end

for i = 1:length(fileList_struct.rawShort)
    fileList_struct.scans{i} = 'all'; %this is a prime example of something that is being kept despite having apparently no function.
end





function [spec,specParams] = GetRawProfileFS_v2(fileList,options,file,DISPLAY_HDRS,X_ZFILLS)
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
spec = [];
specParams = [];

if ~isfield(options,'getSpec'), options.getSpec = 1; end
if ~isfield(options,'sumScans'), error('sumScans not specified'); end

%get Raw and Detector handles
%%%% hRaw = GetRawHandle(fileList.rawFull{file});
%%%% hDetector = GetDetectorHandle(hRaw);

if isunix 
	['wine C:\\python27\\python.exe ', which('MSFileReaderPy.py'), ' -i "' fileList.rawFull{file}, '" -o "', fullfile(fileList.datDir, 'temp.mat"')]
    system(['wine C:\\python27\\python.exe ', which('MSFileReaderPy.py'), ' -i "' fileList.rawFull{file}, '" -o "', fullfile(fileList.datDir, 'temp.mat"')]);
else
    system([which('MSFileReaderPy.exe'), ' -i "' fileList.rawFull{file}, '" -o "', fullfile(fileList.datDir, 'temp.mat"')]);
end

matPy = load(fullfile(fileList.datDir, 'temp.mat'));

%get the list of scans for this spectrum
if isequal(fileList.scans{file},'all')
    %%%% scanList = hDetector.get('FirstSpectrum'):hDetector.get('LastSpectrum');
    scanList = matPy.temp.ScanList;
else
    scanList = fileList.scans{file};
end

%get handle to spectrum and header for each scan
%%%% hFilters = hDetector.get('Filters');
hFilters = matPy.temp.Filters;


%%%% filtersCount = hFilters.count;
filtersCount = matPy.temp.FiltersCount;

%if there is only one filter: easy
if filtersCount == 1
    filter = 1;
    %%%% hFilter = hFilters.Item(filter);
    hFilter = hFilters(1);
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
            hFilter = hFilters(filterIndex);
            %%%% if isequal(hFilter.Text,targetFilter), found=1; break; end
            if isequal(char(hFilter),targetFilter), found=1; break; end
        end
        if ~found, error('Cannot find filter in GetSpectraHandle.m'); end
        filterOfScanIndex(i) = filterIndex;
% THE FOLLOWING TAKES UP TOO MUCH MEMORY HERE
%         hSpectra = hDetector.get('Spectra', filterIndex);
%         hSpectrum(i) = hSpectra.ScanNumber(scanList(i));
%         hSpectrumHeader(i) = hSpectrum(i).get('Header');  %get header for parameters
    end
end

%the list and count of filters actually used is now updated
filtersUsed = unique(filterOfScanIndex);
filtersCount = length(filtersUsed);

%get the centroid (labeldata)
% if options.sumScans
%     %concat the centroids from each scan
%     spectrum.centroidsMZ = [];
%     spectrum.centroidsData = [];
%     for i=1:numberOfScans
%         temp = hSpectrum(i).get('LabelData');
%         spectrum.centroidsMZ = [spectrum.centroidsMZ temp(1,:)];
%         spectrum.centroidsData = [spectrum.centroidsData temp(2,:)];
%     end
%     [spectrum.centroidsMZ,idx] = sort(spectrum.centroidsMZ);
%     spectrum.centroidsData = spectrum.centroidsData(idx);
% elseif options.stitch
    spectrum.centroidsMZ = [];
    spectrum.centroidsData = [];
% else
%     for i=1:numberOfScans
%         temp = hSpectrum(i).get('LabelData');
%         spectrum(i).centroidsMZ = temp(1,:);
%         spectrum(i).centroidsData = temp(2,:);
%     end
% end

%get the spectra
if DISPLAY_HDRS, disp(['File: ',fileList.rawShort{file},'...']); end
tic; lastT=0;

%loop through each filter (these will form the spectra)
spectraScans = {};
spectraFilter = {};

%options.sumScans = 1;
if options.sumScans
    
    for k = 1:filtersCount

        filter = filtersUsed(k);
        scansThisFilter = [];
        d = [];

        %first get a list of scans from this filter
        idx = find(filterOfScanIndex==filtersUsed(k));
        filterScans = scanList(idx);
        scanIT = [];
        scanTIC = [];

        %get the handle to the spectra for this filter
        %%%% hSpectra = hDetector.get('Spectra', filter);

        %loop through each scan in this filter
        for i=1:length(filterScans)

            %get the spectrum for this scan
            %%%% hSpectrum = hSpectra.ScanNumber(filterScans(i));
            %%%% hSpectrumHeader = hSpectrum.get('Header');

            %first check for varying conversion parameters by extracting the trailer extra info
            %%%% hTrailerExtras = hDetector.TrailerExtras();
            %%%% hTrailerExtra = hTrailerExtras.Item(filterScans(i));
            %%%% trailerExtraData = hTrailerExtra.get('Data');
            %%%% newParams.A = str2num(trailerExtraData{2,25})*1e3;
            newParams.A = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.A'))*1e3;
            if i==1, fParams.A = newParams.A;
            elseif newParams.A~=fParams.A
                warning(['Parameter A is changing in ',fileList.rawShort{file},' - included, but check scan ',num2str(filterScans(i))]);
            end
            
            %%%% newParams.B = str2num(trailerExtraData{2,26})*1e6;
            newParams.B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
            if i==1, fParams.B = newParams.B;
            elseif newParams.B~=fParams.B
                warning(['Parameter B is changing in ',fileList.rawShort{file},' - included, but check scan ',num2str(filterScans(i))]);
            end

            
            %%%% scanTIC = [scanTIC hSpectrumHeader.TIC];
            scanTIC = [scanTIC eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.TIC'))];

            %%%% scanIT = [scanIT str2num(trailerExtraData{2,3})];
            scanIT = [scanIT eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'))];
            % eval(strcat('matPy.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'));
            
            %%%% if str2num(trailerExtraData{2,3}) >= 750
            if eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT')) >= 750
                warning(['Underfill in ',fileList.rawShort{file},' - NOT!!! ignoring scan ',num2str(filterScans(i))]);
                %continue;
            end
            
            if options.getSpec
                %get spectrum
                %%%% tempSpectrum = hSpectrum.get('Data');

                %can store the data values as singles to save space - can't do the same with m/z values, need the precision

                %%%% scanmz = tempSpectrum(1,:);
                scanmz = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scanmz'));
                %%%% scand = tempSpectrum(2,:);
                scand = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scand'));
                
                % to extract the noise data, one would used
                % hSpectrum.get('NoiseData');
                % which returns 4 rows: m/z, noise and baseline.
                
                %if the current A or B parameters have changed (are not the selected global parameters), normalise the m/z values to the global A and B parameter
                if newParams.B~=fParams.B || newParams.A~=fParams.A
                    scanmz = mz2f(scanmz,newParams);
                    scanmz = f2mz(scanmz,fParams);
                end
                
                %store scan data
                if isempty(d)
                    d = scand;
                    mz = scanmz;
                else    %append data
                    newd = [];
                    for j=1:size(d,1)
                        %in the case that there are multiple scans, we need to keep the
                        %number of data point entries the same, such that all d's can
                        %use the same 'mz'
                        [newd(j,:),scand_temp,newmz] = FillOut(d(j,:),mz,scand,scanmz);
                    end
                    d = sum([newd;scand_temp],1);   %sum here to increase efficiency
                    mz = newmz;
                end
            end

            scansThisFilter = [scansThisFilter filterScans(i)];

        end

        spectraScans{k} = scansThisFilter;
        spectraFilter{k} = filter;

        %average the parameters if summing spectra over each filter: save each value otherwise
        IT(k) = mean(scanIT);
        TIC(k) = mean(scanTIC);
        
        %%%% mzStart(k) = hSpectrumHeader.LowMass;
        mzStart(k) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.LowMass'));
		
        %%%% mzEnd(k) = hSpectrumHeader.HighMass;
        mzEnd(k) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.HighMass'));
        
        params(k) = fParams;
%         else
%             IT((end+1):specCount) = scanIT;
%             TIC((end+1):specCount) = scanTIC;
%             mzStart((end+1):specCount) = hSpectrumHeader.LowMass;
%             mzEnd((end+1):specCount) = hSpectrumHeader.HighMass;
%             params((end+1):specCount).A = fParams.A;
%             params((end+1):specCount).B = fParams.B;
%         end

        if options.getSpec
            %average spectra abundance
            data = d ./ length(scansThisFilter);
            %convert to frequency (increasing values); decalibrate
            spectrum(k).f = mz2f(mz,params(k));
            spectrum(k).f = fliplr(spectrum(k).f);
            spectrum(k).data = fliplr(data);
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

else    %not summing scans

    scanCount = 0;;
    for k = 1:filtersCount

        filter = filtersUsed(k);

        %first get a list of scans from this filter
        idx = find(filterOfScanIndex==filtersUsed(k));
        filterScans = scanList(idx);

        %get the handle to the spectra for this filter
        %%%% hSpectra = hDetector.get('Spectra', filter);

        %loop through each scan in this filter
        for i=1:length(filterScans)

            %get the spectrum for this scan
            %%%% hSpectrum = hSpectra.ScanNumber(filterScans(i));
            %%%% hSpectrumHeader = hSpectrum.get('Header');
            
            scanCount = scanCount + 1;
            spectraScans{scanCount} = filterScans(i);
            spectraFilter{scanCount} = filter;

            %%%% hTrailerExtras = hDetector.TrailerExtras();
            %%%% hTrailerExtra = hTrailerExtras.Item(filterScans(i));
            %%%% trailerExtraData = hTrailerExtra.get('Data');
            
            %%%% TIC(scanCount) = hSpectrumHeader.TIC;
            TIC(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.TIC'));
            %%%% IT(scanCount) = str2num(trailerExtraData{2,3});
            IT(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'));
            %%%% mzStart(scanCount) = hSpectrumHeader.LowMass;
            mzStart(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.LowMass'));
            %%%% mzEnd(scanCount) = hSpectrumHeader.HighMass;
            mzEnd(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.HighMass'));
            %%%% params(scanCount).A = str2num(trailerExtraData{2,25})*1e3;
            params(scanCount).A = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.A'))*1e3;
            %%%% params(scanCount).B = str2num(trailerExtraData{2,26})*1e6;
            params(scanCount).B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
            
            if options.getSpec
                %get spectrum
                %%%% tempSpectrum = hSpectrum.get('Data');

                %can store the data values as singles to save space - can't do the same with m/z values, need the precision
                %%%% scanmz = tempSpectrum(1,:);
                scanmz = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scanmz'));
                %%%% scand = tempSpectrum(2,:);
                scand = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scand'));

                %progress counter
                progress = i/length(filterScans)*100;
                t=toc;
                if t<5, lastT = floor(t/10);
                elseif t>5 && floor(t/10)>lastT
                    disp([num2str(progress),'%']);lastT=floor(t/10);drawnow;
                end

                %convert to frequency (increasing values); decalibrate and store
                spectrum(i).f = mz2f(scanmz,params(i));
                spectrum(i).f = fliplr(spectrum(i).f);
                spectrum(i).data = fliplr(scand);
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
        T(i) = 1/deltaf/2^X_ZFILLS;
        if round(T(i)*1000) ~= 768, warning(['T not 0.768 in raw file, is ',num2str(T(i))]); keyboard; end;

        %each data point will have a new index relative to the first frequency
        newIdx = (spectrum(i).f-spectrum(i).f(1))./deltaf;  %index relative to first freq. val
        if find(abs(round(newIdx)-newIdx)>0.3), warning(['Possible error finding data indices in GetRawProfileFS.m! Press any key to continue...']); end %pause; end
        newIdx = round(newIdx) + 1;

        %spread out the data into the new data matrix, at the same time padding the gaps with zero
        spectrum(i).newData(newIdx) = spectrum(i).data;

        %the new frequency values are the indices multiplied by the frequency step, plus the starting frequency
        spectrum(i).newf = ((0:length(spectrum(i).newData)-1).*deltaf) + spectrum(i).f(1);
    end
end

%return values of interest
numberOfSpectra;
for i=1:numberOfSpectra
    if options.getSpec
        spec(i).y = spectrum(i).newData;
        spec(i).f = spectrum(i).newf;
        specParams(i).T = T(i);
    end
%     spec(i).xCentroidsMZ = spectrum(i).centroidsMZ;
%     spec(i).xCentroidsData = spectrum(i).centroidsData;
    specParams(i).mzStart = mzStart(i);
    specParams(i).mzEnd = mzEnd(i);
    specParams(i).TIC = TIC(i);
    specParams(i).A = params(i).A;
    specParams(i).B = params(i).B;
    specParams(i).C = 0;
    specParams(i).zFills = X_ZFILLS;
    specParams(i).IT = IT(i);
    specParams(i).filterIndex = spectraFilter{i};
    specParams(i).scans = spectraScans{i};
end

% release all interfaces derived from the Raw COM server
% and then delete the server itself.
%%%% hRaw.delete;


function [transient,specParams,message] = ReadInDat_3_0(fDir,fStem,scanList,n,PARAMS,DISPLAY_HDRS)
%read in n data points from a list of .dat transient files in scanList
% 
%given:
% fStem file stem
% n, number of data points to read in (n = 'all' -> read all data available)
% scanlist: if scanList is empty, read in all available transients
% if PARAMS.IGNORECHANGINGCALPARAMS_ON=1, transients with different
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
            if PARAMS.IGNORECHANGINGCALPARAMS_ON
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

