function SampleFilter(pls_dso_file, filt_thold, class_name, dso_outfile)
%
%
%   SampleFilter_BB(pls_dso_file,filt_thold,class_name,dso_outfile)
%
%   This function filters a PLS Toolbox dataset object, applying a
%   threshold (filt_thold e.g. 80) to the peaks within a 'class' of data.
%
%   Any peak that is above the threshold in ANY class will be retained in
%   the output dataset.
%
%   In order to use this script with non-representative groups e.g. pooled
%   QCs and blanks, the script will ignore any spectra/classes that are not
%   checked as 'include'. 
%
%   Suggested use - run 'AlignSample' then 'CreateDSO' which will
%   produce a large, unfiltered dataset object. Save this to disk somewhere, then open in
%   Matlab by 'double-clicking' the dataset in the workspace, label classes, deselect any QCs etc, save
%   in the workspace, close the dataset GUI. Then run this script to
%   filter.
%
%   Inputs: 
%       fileList_in:	Full path to an xml representation of a fileListManager object containing various common variables.
%       pls_dso_file: 	Full path to an xml representation of a PLS Toolbox dataset object.
%	filt_thold:	the filter threshold. e.g. 80
%	class_name:	Either a string representing the label for a list of classes (within the dataset object) or an empty string '' if peaks are to be filtered across ALL classes.
%	dso_outfile:	Full path to an xml file for the xml dataset object to be written to
%
%
%   Outputs:
%	this tool returns no Matlab variables, only writing to  the file specified by dso_outfile.
%
%   R.L.Davidson 06/11/2014
%
%   Version 3.3

%% CHECK FOR PLS TOOLBOX IN PATH


if isa(filt_thold, 'char')
	filt_thold = str2num(filt_thold);
end

if ~isdeployed
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


%% PARAMS


% check for use of class specific filtering

class_name_flag = 1;

if isempty(class_name) %check for an empty class name  
	class_name_flag = 0;
end




%% IMPORT PLS Toolbox Dataset Object
[dso, name, source] = autoimport(pls_dso_file, 'xml');
%% VALIDATE DATA AND PARAMS

if class_name_flag
    
    class_list_index = [];
    for i = 1:size(dso.classname,2)
        if(strcmp(dso.classname{1,i},class_name))
            class_list_index(end+1) = i;
        end
    end
    
    if length(class_list_index)>1
        sprintf('error: more than one class with the name %s',class_name)
        dso_out = [];
        unique_peaks = [];
        return
    elseif isempty(class_list_index)
        sprintf('error: no class with the name %s',class_name)
        dso_out = [];
        unique_peaks = [];
        return
    end
    
    

    classes = dso.class{1,class_list_index};
    class_lookup = dso.classlookup{1,class_list_index};

else
    classes = ones(1,size(dso.data,1));
    class_lookup = {1 'all_samples'};
end

%% EXTRACT DATA FROM DSO

data = dso.data;
includes = dso.include{1};
mz = dso.axisscale{2};
if isempty(mz)
    mz = 1:size(dso.data,2);
end
name = dso.name;
labels = dso.label{1};
class_IDs = unique(classes(includes));

%clear dso  %retain for filling in classlabels in dso_out etc later. 

%% FILTER EACH CLASS

for i = 1:length(class_IDs)
    tmp_ind = includes(find(classes(includes)==class_IDs(i)));
    tmp_data = data(tmp_ind,:);
    tmp_data(tmp_data>0)=1;
    tmp_data_sums = sum(tmp_data,1);
    
    tmp_thold = ceil(size(tmp_data,1)/100*filt_thold);
    
    keep_ind{i} = find(tmp_data_sums>=tmp_thold);
    totals{i} = tmp_data_sums;
    percents{i} = tmp_data_sums./size(tmp_data,1);
    
end

clear tmp_ind tmp_data tmp_data_sums 

total_keep_ind = [];
for i = 1:length(keep_ind)
    total_keep_ind = [total_keep_ind,keep_ind{i}];
end
total_keep_ind = unique(total_keep_ind);

%% REDUCE DATA

data = data(:,total_keep_ind);
mz = mz(total_keep_ind);

dso_out = dataset(data);
unique_peaks.mz_order = mz';
%% PICK OUT UNIQUE PEAKS

for i = 1:length(class_IDs)
    unique_peaks.class_order{i} = class_lookup{find([class_lookup{:,1}]==class_IDs(i)),2};
    
    unique_peaks.num_peaks(:,i) = totals{i}(total_keep_ind);
    unique_peaks.percent(:,i) = percents{i}(total_keep_ind);
end

% from original version of samplefilter_3_0
%{
%map each 'keep index' onto the new total_keep_ind
for i = 1:length(keep_ind)
   for j = 1:length(keep_ind{i})
       keep_ind{i}(j) = find(keep_ind{i}(j)==total_keep_ind);
   end
end

for i = 1:length(keep_ind)
    tmp_val = keep_ind{i};
    tmp_range = 1:length(keep_ind);
    tmp_range(i) = [];
    tmp_val2 = [];
    for j = tmp_range
        very_tmp_ind = keep_ind{j};
        tmp_val2 = [tmp_val2,very_tmp_ind];
    end
    tmp_val2 = unique(tmp_val2);
    tmp_unique_ind = setdiff(tmp_val, tmp_val2);
    unique_peaks(i).class = class_lookup(find([class_lookup{:,1}]==class_IDs(i)),2);
    unique_peaks(i).mz = mz(tmp_unique_ind);
    unique_peaks(i).ind = tmp_unique_ind;

end
%}


%% COPY META INFO FROM DSO TO DSO_OUT
dso_out.include{1} = includes;
dso_out.axisscale{2} = mz;


for i = 1:size(dso.labelname,1)
    for j = 1:size(dso.labelname,2)
        dso_out.labelname{i,j} = dso.labelname{i,j};
    end
end

for i = 1:size(dso.classname,1)
    for j = 1:size(dso.classname,2)
        dso_out.classname{i,j} = dso.classname{i,j};
    end
end

for i = 1:size(dso.label,1)
    for j = 1:size(dso.label,2)
        dso_out.label{i,j} = dso.label{i,j};
    end
end

for i = 1:size(dso.class,1)
    for j = 1:size(dso.class,2)
        dso_out.class{i,j} = dso.class{i,j};
    end
end

for i = 1:size(dso.classlookup,1)
    for j = 1:size(dso.classlookup,2)
        dso_out.classlookup{i,j} = dso.classlookup{i,j};
    end
end

dso_out.name = [dso.name,'_samplefiltered',num2str(filt_thold)];



%export DSO to XML
autoexport(dso_out, dso_outfile, 'xml');



%% RETURN PATH TO ORIGINAL STATE
%clear classes (required, but only to be implemented after file is saved somewhere...)
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end





