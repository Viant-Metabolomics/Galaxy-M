function CreateDSO_BB(fileList_in, matrix_file, row_file, col_file, outfile)
%
%   CreateDSO(fileList_in, matrix_file, row_file, col_file, outfile)
%
%   This function loads the files in the 'combined peaks' folder created by
%   AlignSample(), and converts the MZ list and Intensity matrix to a PLS
%   Toolbox dataset object.
%
%   NB: The dataset object is NOT saved anywhere! You must do this from
%   within matlab.
%
%   Inputs:
%       1, fileList_in	:the full path to an xml fileList 
%       2, matrix_file	:the full path to a peak matrix .txt file as produced by AlignSample
%	3, row_file	:the full path to a row headings file (i.e. list of sample names) as produced by AlignSample
%	4, col_file	:the full path to a column headings file (i.e. list of MZ values) as produced by AlignSample
%	5, outfile	:the full path to an xml file for storing the xml representation of the created PLS Toolbox dataset object
%
%   Outputs:
%       No Matlab variables are returned by this tool, the only output being the xml file written to 'outfile'.
%
%
%   R.L.Davidson 26/11/12
%
%   Version 3.0

if ~isdeployed
	%% CHECK FOR PLS TOOLBOX IN PATH
	try
		dtst = dataset(rand(10,100));       %attempt to create a dataset object
		props = properties(dtst);           %request the properties of said object
	catch err
		if strcmp(err.identifier,'MATLAB:UndefinedFunction') %if dataset function not available then neither PLS Toolbox or Stats Toolbox are installed.
		    sprintf('Matlab does not recognise dataset function. Neither PLSToolbox or Stats toolbox installed. Please amend and try again.')
		    return
		end
	end


	if ~isempty(props) %PLS Toolbox datasets have no properties at initiation whereas Matlab datasets have 2.
		
		% if here, Statistics Toolbox version has been used. need to move stats toolbox below pls toolbox.
		clear dtst
		clear classes   %need to remove the statistics toolbox dataset class
		
		original_path = path; %save original path
		rem = original_path;
		pls_path = '';
		rem_path = '';
		while true
		    [str,rem] = strtok(rem,pathsep);
		    if isempty(str)
		        break
		    elseif strfind(str,'pls_toolbox') %covers all pls_toolbox entries
		        pls_path = [pls_path,str,pathsep];
		    else
		        rem_path = [rem_path,str,pathsep];
		    end
		end
		
		if ~isempty(pls_path) %check for no PLS toolbox installed
		    path(pls_path,rem_path); %put PLS at the top!
		    rehash pathreset;
		    rehash toolboxreset;
		    
		    dtst = dataset(rand(10,100));
		    props = properties(dtst);
		    if ~isempty(props) % if rehash has not worked, quit.
		        sprintf('Cannot appropriately rejig path. Please manually place PLSToolbox above Stats Toolbox in path.')
		        path(original_path);
		        return
		    end
		    
		else    % If no stats toolbox entries have been found in path there is a more serious problem.
		    sprintf('PLS Toolbox not on path. Please Install and try again.')
		    return
		end
		
	else
		sprintf('PLS Toolbox dataset objects are available. Continuing.')
		original_path = path;
	end
end

%% PARAMETERS
BLANKSTR = 'blank';

%% LOAD DATA FROM FILELIST

fileList = import_filelist_xml(fileList_in);

PARAMS.NUMREPS = fileList.numReps;     % the default number of replicates

%% LOAD DATA FROM txt FILES


try mz = load(col_file);
catch
	sprintf('ERROR: could not load col headers: %s',mz_file)
        return
end
    

try peakMatrix = load(matrix_file);
catch
	sprintf('ERROR: could not load intensities: %s',int_file)
        return
end
    

try fid = fopen(row_file);
catch
        sprintf('ERROR: could not load row headers: %s',row_file)
        return
end
    
count=1;
while 1
        tline = fgetl(fid);
        if ~ischar(tline),break,end
        sample_name{count,1} = tline;
        count = count+1;
end
fclose(fid);

%count = 0;
classnames_BF = {};
for si=1:PARAMS.NUMREPS:length(fileList.rawShort)
    %count = count+1;
    %sampleGpPks.blankflag(count) = strncmpi(fileList.sampleID{si},BLANKSTR,length(BLANKSTR));
    if strncmpi(fileList.spec(si).sampleID,BLANKSTR,length(BLANKSTR));
        classnames_BF{end+1} = 'blank';
    else
        classnames_BF{end+1} = 'sample';
    end
end

classnames = {};
for si=1:PARAMS.NUMREPS:length(fileList.rawShort)
    classnames{end+1} = fileList.spec(si).sampleID;
end

dso = dataset(peakMatrix);
dso.axisscale{2} = mz;

for i = 1:length(sample_name)
    ind = strfind(sample_name{i},fileList.setUID);
    ind = ind+length(fileList.setUID)+1; %removes the dataset name and the proceeding '_'
    sample_name{i} = sample_name{i}(ind:end-4); %end-4 removes the .txt at the end
end

dso.label{1} = sample_name;
dso.labelname{1} = 'sample_name';
%{
dso.class{1} = sampleGpPks.blankflag+0;
dso.classname{1} = 'blank_flag';
dso.classlookup{1,1}.assignstr = {0 'sample'};
dso.classlookup{1,1}.assignstr = {1 'blank'};
%}
dso.classname{1,1} = 'classes';
dso.classid{1,1} = classnames;

dso.classname{1,2} = 'blank_flag';
dso.classid{1,2} = classnames_BF;

dso.name = [fileList.setUID,'_dso'];

%% Save DSO to DSO Directory

%DSODir = fileList.DSODir;
%mkdir(DSODir); %should be the first instance of the directory being created...
%DSO_save_path = [DSODir,'dso_Create.mat'];
%save(DSO_save_path,'dso');

%make an output for Galaxy
%fid = fopen(outfile,'w');
%fprintf(fid,'CreateDSO - Success!');
%fclose(fid);

%export DSO to XML
autoexport(dso, outfile, 'xml');

%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end



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
    sprintf('Failed to read XML file.','Import error');
end

if ~tree.hasChildNodes()
    sprintf('XML tree empty!','Import error');
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
    sprintf('XML Document root is not fileList - stopping')
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
    sprintf('XML file does not contain all fields - stopping!','Import error')
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
    sprintf('XML fileList does not contain information about spectral raw files... quitting','Import error')
    return
end

%add some conversions to help this xml format fit with the expected
%fileList structure variable names from legacy code:

%SumTransients conversions
fileList_struct.avTransDir = fileList_struct.avgdTransDirectory;
fileList_struct.setUID = fileList_struct.setIdentifier;

fileList_struct.rawShort = {fileList_struct.spec.ID};
for i=1:length(fileList_struct.rawShort)
    fileList_struct.rawFull{i} = [fileList_struct.rootDirectory,fileList_struct.rawDirectory, fileList_struct.rawShort{i},'.raw'];
end
fileList_struct.rawContig = 0; % set to 1 if raw files all combine to form a single spectrum (never used)

fileList_struct.datDir = [fileList_struct.rootDirectory, fileList_struct.datDirectory];
for i=1:length(fileList_struct.rawShort)
    fileList_struct.datStem{i} = [fileList_struct.rawShort{i},'_'];
end

for i = 1:length(fileList_struct.rawShort)
    fileList_struct.scans{i} = 'all'; %this is a prime example of something that is being kept despite having apparently no function.
end

return
end



