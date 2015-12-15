function CreateDSO(fileList_in, matrix_file, row_file, col_file, outfile)
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

fileList = ImportFileListXML(fileList_in);

PARAMS.NUMREPS = fileList.nReplicates;     % the default number of replicates

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
classnames_Batch = {};
classnames_Order = {};
for si=1:PARAMS.NUMREPS:length(fileList.Samples)
    if strncmpi(fileList.Samples(si).sampleID,BLANKSTR,length(BLANKSTR));
        classnames_BF{end+1} = 'blank';
    else
        classnames_BF{end+1} = 'sample';
    end
    classnames_Batch{end+1} = fileList.Samples(si).batchID;
    classnames_Order{end+1} = fileList.Samples(si).orderID;
end

classnames = {};
for si=1:PARAMS.NUMREPS:length(fileList.Samples)
    classnames{end+1} = fileList.Samples(si).sampleID;
end

dso = dataset(peakMatrix);
dso.axisscale{2} = mz;

for i = 1:length(sample_name)
    %ind = strfind(sample_name{i},fileList.setUID);
    %ind = ind+length(fileList.setUID)+1; %removes the dataset name and the proceeding '_'
    sample_name{i} = sample_name{i}(7:end-4);
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

dso.classname{1,3} = 'batch_info';
dso.classid{1,3} = classnames_Batch;

dso.classname{1,4} = 'order_info'; % Order of measurements
dso.classid{1,4} = classnames_Order;

dso.name = ['dso'];

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



