function BlankFilter_BB(pls_dso_file, thresh_perc, maxmean, thresh_fold, outfile)
%
% BlankFilter_BB( pls_dso_file,thresh_perc ,maxmean , thresh_fold, output)
%
%   inputs- a .csv file with the following parameters:
%	
%	.pls_dso_file -	    Full path to an xml representation of the PLS Toolbox dataset object as created by CreateDSO. This contains %				certain class labels etc that are important for these filters.	
%       .thresh_perc -      The percentage of samples that must beat the blank
%                               filter. This test only compares those samples that
%                               have an intensity. 
%                               So, if you choose 50, for each peak that is in a
%                               blank, 50% of the samples that have a matching peak must have a 'better'
%                               value for the whole peak column to be kept.
%                               Otherwise all peaks are removed! DEFAULT, 100
%       .maxmean -          Set to 0 to use the max of the blanks, 1 to use the mean
%                               When using the mean, missing values are counted as
%                               zero. (this is how the old script used to
%                               do it) DEFAULT, 1 (mean)
%       .thresh_fold -      How much bigger do you want the samples to be when compared to the blank?                       
%                               For samples to be 2 and a half times larger than
%                               blanks, set to 2.5 - for 10 times larger,
%                               set to 10. DEFAULT, 10
%	.outfile - 	    Full path to an xml file for receiving an xml representation of the filtered PLS Toolbox dataset object.
%
%   outputs:
%       No matlab variables are returned by this function, the only output being the dataset object in xml format that is written to 'outfile'
%
%   R.L.Davidson 08/02/2013
%
%   Version 3.0


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



%% LOAD DSO


[dso, name, source] = autoimport(pls_dso_file, 'xml');


%% VALIDATE DATA AND PARAMS
blank_class_index = 0;
for i = 1:size(dso.classname,2)
if strcmp(dso.classname{1,i},'blank_flag')
blank_class_index = i;
end
end
yes_blanks = 1; %do some error checking, looking for samples labeled 'blank'
if ~blank_class_index
disp('error: no class list named ''blank_flag''. Please label class list and retry.')
yes_blanks = 0;
end

blank_lookup_index = find(strcmp(dso.classlookup{1,blank_class_index}(:,2),'blank'));

if isempty(blank_lookup_index)
disp('error: blank_flag classes do not contain a class named ''blank''. Please check. ')
yes_blanks = 0;
end

if yes_blanks
    blank_class = dso.classlookup{1,blank_class_index}{blank_lookup_index,1};

    blank_index_s = find(dso.class{1,blank_class_index}==blank_class);

    if isempty(blank_index_s)
        disp('error: no samples are classed as ''blank'' in class ''blank_flag''. No filtering possible. Use dso_in instead of creating dso_out. ')
        yes_blanks = 0;
    end
    
end
yes_samples = 1; %a flag for noting whether samples (non blanks of any label) are 'included'
if yes_blanks
    
    data = dso.data;
    include_peak = dso.include{2};
    include_sample = dso.include{1};
    include_blank_index = find(dso.class{1,blank_class_index}(include_sample)==blank_class);
    include_sample_index = find(dso.class{1,blank_class_index}(include_sample)~=blank_class);
    
    if isempty(include_blank_index)
        disp('error: your blanks are not included. Please open DSO and ''include'' some blanks.')
        yes_blanks = 0;
    elseif isempty(include_sample_index)
        disp('error: you have not included any non-blanks! Please open DSO and ''include'' some samples.')
        yes_samples = 0;
    end
    
end

%% PERFORM FILTER
if yes_blanks & yes_samples
    
data_include = data(include_sample, include_peak);

if maxmean
    blank_comb = mean(data_include(include_blank_index,:),1);
else
    blank_comb = max(data_include(include_blank_index,:),[],1);
end

blank_index_pk = find(blank_comb);

%n_included = size(data_include,1);
%n_thresh = ceil(thresh_perc*(n_included/100));

keep_index_bin = ones(1,size(data_include,2));
for i = 1:size(blank_index_pk,2)
    pks = data_include(include_sample_index,blank_index_pk(i));
    pks = pks(find(pks));
    tmp_thresh = ceil(thresh_perc*(length(pks)/100));
    pass_index = find(  (pks./blank_comb(blank_index_pk(i)))>=thresh_fold);
    if length(pass_index)<tmp_thresh
        keep_index_bin(blank_index_pk(i))=0;
    end
    
end

keep_index_actual = find(keep_index_bin);

%at this point, we've compared the peaks to the blanks... now we need to
%take those results and match them to the original data. 

%keep_index (bin and actual) are matched to data_include. Need to match to
%all data. 

not_include_peak = setdiff(1:size(data,2),include_peak);
keep_include_peak = include_peak(keep_index_actual);

[all_output_peaks,all_output_peak_index] = unique([not_include_peak, keep_include_peak]); %this combines any exluded peaks, and the retained, filtered peaks and sorts the indexes


%% CREATE OUTPUT DATASET

%NB: we will output a DSO with all samples from the original and this
%combined set of included+filtered and not-included peaks. We will set the
%include variables in the DSO to reflect any not-included options from the
%original. 

data_out = data(:,all_output_peaks);
include_peaks_out = 1:size(data_out,2);

%using the sort index from the unique command above, we can check which of
%the output peaks were originally in not_include_peak (i.e. should be
%excluded here)
if ~isempty(not_include_peak)       %check if needed
    exclude_peaks_out = find(all_output_peak_index<=length(not_include_peak)); %find the resorted index positions of peaks from not_include_peak
    include_peaks_out(exclude_peaks_out) = []; %remove those unincluded peaks from this include list!
end

dso_out = dataset(data_out);

dso_out.axisscale{2} = dso.axisscale{2}(all_output_peaks);
dso_out.axisscale{1}  = dso.axisscale{1};

dso_out.include{1} = include_sample;
dso_out.include{2} = include_peaks_out;

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

dso_out.name = [dso.name,'_blankfiltered'];
else
    dso_out = dso;
    dso_out.name = [dso.name, '_no_blankfilter'];

end


%% SAVE OUTPUT DATASET
%dso = dso_out;
%DSO_save_path = [DSODir,'dso_BlankFilter.mat'];
%save(DSO_save_path,'dso');

%make an output for Galaxy
%fid = fopen(outfile,'w');
%fprintf(fid,'BlankFilter - Success!');
%fclose(fid);

%export DSO to XML
autoexport(dso_out, outfile, 'xml');

%% RETURN PATH TO ORIGINAL STATE
%clear classes (required, but only to be implemented after file is saved somewhere...)
path(original_path);
rehash pathreset;
rehash toolbox;

return
end



