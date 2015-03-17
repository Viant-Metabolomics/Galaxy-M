function GetTopPCALoadings(file, lv_choice, top_pcent,results_file) 
%
% GetTopPCALoadings(file, lv_choice, top_pcent,results_file)
%
%   This function retrieves the top X percent of variables (metabolite
%   peaks) based on the absolute values of the PCA loadings of a chosen
%   Latent Variable. 
%
%   It returns the M/Z value and average intensity of each
%   variable in a file that is intended for submission to the MI Pack
%   metabolite identification package.
%
%
%   Inputs:
%       1, file	:       the full path to an xml representation of a PLS Toolbox model (PCA) 
%       2, lv_choice:   the desired latent variable (integer)
%       3, top_pcent:   a percentage value for thresholding the loadings. 
%       4, results_file:the full path to a tab delimited file that will
%       receive the output (a two column list of M/Z and average Intensity
%       values, with headings M/Z and INTENSITY)
%       
%   Outputs:
%       No Matlab variables are returned by this tool, the only output being the tab delimited file written to 'results_file'.
%
%
%   R.L.Davidson 12/12/2014
%
%   Version 1.0


%% BASIC ERROR CHECKING

if top_pcent>100
    error('Maximum threshold is 100 percent of variables, %d chosen.',ceil(top_pcent))
end

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



%% IMPORT MODEL

pca_model = autoimport(file, 'xml');

loadings_dim = size(pca_model.loads{2});
if lv_choice>loadings_dim(2)
    error('PCA Model does not have more than %d LVs: requested LV was %d.',loadings_dim(2),lv_choice);
end

lv = pca_model.loads{2}(:,lv_choice);
mz = pca_model.detail.axisscale{2}(1,:);


include_1 = pca_model.detail.includ{1};
include_2 = pca_model.detail.includ{2};
data = pca_model.detail.data{1}.data(include_1,include_2);

dim = size(data);

if dim(1)>1
    intensity = mean(data);
else
    intensity = data;
end

[lv_sort,lv_origin] = sort(abs(lv),'descend'); 

threshold = ceil(length(lv_sort)/100*top_pcent);

top_index = sort(lv_origin(1:top_pcent), 'ascend');
top_mz = mz(top_index);
top_intensity = intensity(top_index);

results = [top_mz',top_intensity'];

fid = fopen(results_file,'w');
if fid == -1
    error('could not open results output file')
end

fprintf(fid, 'M/Z\tINTENSITY\n');
fprintf(fid, '%f\t%f\n',results');
fclose(fid);


