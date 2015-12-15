function [dtst_out,time_taken] = ImputeKNN(dso_xml, k, col_t, row_t, outfile)
%
%   [dtst_out,time_taken] = impute_knn(dtst,k, col_t, row_t)
%
%   inputs:
%       dtst = the matrix or dataset object (samples in rows)
%       k = how many neighbours to use to estimate missing values
%       col_t = the min percentage of rows that must contain each peak e.g. 50 
%       row_t = the max percentage of missing values any one row can have
%       e.g. 30
%
%   outputs:
%       dtst_out = the matrix or dataset object with missing values
%       estimated
%       time_taken = the time in minutes that the calculation has taken
%       (only for testing purposes)
%
%
%   R.L.Davidson 10/01/2012
%
%*****************************************************************************

if isa(k, 'char')
	k = str2num(k);
end

if isa(col_t, 'char')
	col_t = str2num(col_t);
end

if isa(row_t, 'char')
	row_t = str2num(row_t);
end


%% CHECK FOR PLS TOOLBOX IN PATH
if ~isdeployed
	try
	    dtst = dataset(rand(10,100));       %attempt to create a dataset object
	    props = properties(dtst);           %request the properties of said object
	catch err
	    if strcmp(err.identifier,'MATLAB:UndefinedFunction') %if dataset function not available then neither PLS Toolbox or Stats Toolbox are installed.
		
		%fid = fopen(messagefile,'w');
		sprintf('Matlab does not recognise dataset function. Neither PLSToolbox or Stats toolbox installed. Please amend and try again.');
		%fclose(fid);
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
		    %fid = fopen(messagefile,'w');
		    sprintf('Cannot appropriately rejig path. Please manually place PLSToolbox above Stats Toolbox in path.');
		    %fclose(fid);
		    path(original_path);
		    return
		end
		
	    else    % If no stats toolbox entries have been found in path there is a more serious problem.
		%fid = fopen(messagefilefile,'w');
		sprintf('PLS Toolbox not on path. Please Install and try again.');
		%fclose(fid);
		path(original_path);
		return
	    end
	    
	else
	    sprintf('PLS Toolbox dataset objects are available. Continuing.')
	    original_path = path;
	end
end



%% LOAD DATASETOBJECT

[dtst, name, source] = autoimport(dso_xml, 'xml');


start_time = cputime;
%check for dataset object
if isobject(dtst)
    x = dtst.data;
    %test for those that are 'included'
    %inc = dtst.include{1};
    %x = x(inc,:);
else
    x = dtst;
end

[m,n] = size(x);


%% PERFORM KNN IMPUTATION

% open  'outfile' for status reporting
%fid = fopen(messagefile,'a');


%{
%should test for missing values in each variable in case there are too
%many. This throws a warning. Requires 'colmax' or col_threshold (col_t)
%}
%check for columns(peaks) where there are less than col_t  entries

min_peaks = floor(m/100*col_t);
col_warn = 0;
col_warn_ind = [];
for i = 1:n
    num_peaks_i = length(find(x(:,i)));
    if num_peaks_i<min_peaks
        col_warn = col_warn+1;
        col_warn_ind(end+1) = i;
    end
end

if col_warn >0
    sprintf('%i columns (peaks) did not meet the required threshold.\n',col_warn);
    %x(:,col_warn_ind) = [];
end

%the R script checks the number of missing values in each row and if it
%falls below a threshold, the mean of each variable is used rather than
% looking for similar expressions in other variables. Requires 'rowmax' or
% row_t

max_missing = floor(n/100*row_t);
row_warn = 0;
row_warn_ind = [];

for i = 1:m
    num_missing_i = length(find(x(i,:)==0));
    if num_missing_i>max_missing
        row_warn = row_warn+1;
        row_warn_ind(end+1) = i;
    end
end

if row_warn >0
    sprintf('%i rows (samples) had too many missing values. For each missing peak, the average of other instances of that peak will be used as the imputed value.\n',row_warn);
end

%remove the rows and columns that are too sparse
impute_col_ind = 1:n;
impute_col_ind(col_warn_ind) = [];

impute_row_ind = 1:m;
impute_row_ind(row_warn_ind) = [];

%create 2 matrices for imputation - 1 for knn, one for averaging. 
x_imputing = x(impute_row_ind, impute_col_ind);
x_averaging = x(row_warn_ind,impute_col_ind);


x_ave_flag = 1;
%check to see if either matrix is empty
if isempty(x_averaging)
    sprintf('No peaks will be averaged, all will be knn imputed (good)\n');
    x_ave_flag = 0;
end
x_imp_flag = 1;
if isempty(x_imputing)
    sprintf('No peaks will be knn imputed! (all will be averaged) Too many missing values in each row!(bad)\n')
    x_imp_flag = 0;
end

if x_ave_flag
    % where the row has too few entries, each peak will have any missing
    % values set to the average for that peak. 
    for i = 1:size(x_averaging,2)
        ave = mean(x_averaging(find(x_averaging(:,i)),i));
        x_averaging(find(x_averaging(:,i)==0),i) = ave;
    end
end


if x_imp_flag
    x_imputed = x_imputing;
    for i = 1:size(x_imputing,1)
        for j = 1:size(x_imputing,2)
            if x_imputing(i,j) ==0
                ind_known = find(x_imputing(:,j));
                ind_full = [i;ind_known]; %an index of all known values in this column plus the unknown value.
                
                x_tmp = x_imputing(ind_full,:); %reduce the matrix to just the desired rows.
                
                z_ind = find(min(x_imputing(ind_full,:))); %find all those peaks that have no zero values in the rows in question
             
                
                %check that there are enough peaks (>=k) that have these values                %values
                if isempty(z_ind)
                    sprintf('Missing value %i, %i cannot be imputed because no other peaks share its values.\n',[i,j]);
                    break
                elseif length(z_ind)<k
                    sprintf('Missing value %i, %i cannot be imputed because less than k(=%i) other peaks share its values.\n',[i,j,k]);
                    break
                else
                    dists = [];
                    for z = z_ind
                        %find the distance between peak 'j' and the peaks
                        %in the 'z' index.
                        dists(end+1) = pdist(x_imputing(ind_known,[j,z])');
                    end
                    
                    %sort the dists to find the k nearest neighbours
                    
                    [sorted_dists, dist_index] = sort(dists, 'ascend');
                    
                    total_distance = sum(sorted_dists(1:k));
                    tmp = 0;
                    for iter = 1:k
                        weight_tmp = 1-sorted_dists(iter)/total_distance;
                        tmp = tmp+ weight_tmp*x_imputing(i,z_ind(dist_index(iter)));
                    end
                    x_imputed(i,j) = tmp/(k-1);
                    
                    
                    
                end
            end
        end
    end
end

    


end_time = cputime;
time_taken = (end_time-start_time)/60;

%% CREATE DATASET FOR OUTPUT


dtst_out = zeros(size(x));
dtst_out(impute_row_ind, impute_col_ind) = x_imputed;
dtst_out(row_warn_ind,impute_col_ind) = x_averaging;

if isobject(dtst) %return a dataset object if that was used initially. 
    dtst_out_object = dtst;
    dtst_out_object.data = dtst_out;
    dtst_out = dtst_out_object;
end



%% SAVE OUTPUT DATASET
autoexport(dtst_out, outfile, 'xml');


%close the messagefile output
%fclose(fid);


%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end
                




