function GetAxisscale(dso_xml,dimension, outfile)
% -------------------------------------------------------------------
% GetAxisscale(dso_xml,dimension, outfile)
% 
% This is a function designed to work within the Birmingham Metabolomics Pipeline in Galaxy.  
%
% It extracts the axisscale in the given dimension, intended to be the MZ values of your mass spectra. This is intended for interoperability purposes - allowing the data to be passed to other pipelines that do not use the PLS Toolbox Dataset Object. 
%
% It makes use of the 'include' function of the dataset object. Changing the 'include' parameter for the Dataset will alter the output of this tool. 
%
% The output is a .csv file with one row (for dimension 2 = MZ axis) or one column (for dimension 1 = Samples), and no headers.
%
% The option is given to access the axisscale for Samples. This could be used in a regression problem where the samples were associated with a continuous variable rather than a discrete classification. 
%
%
% inputs:
%	dso_xml	:	the path to an xml format PLS Toolbox Dataset Object as output by autoexport()
%	dimension:	an integer that can be either 1 or 2, representing Sample Labels or MZ respectively.
%	outfile	:	the path to a .csv file for outputting the two columns
% -------------------------------------------------------------------

%% SANITY CHECK DIMENSION
if dimension <1 | dimension >2
    sprintf('Dimension inappropriate should be 2 or 1')
    return
end

%% CHECK FOR PLS TOOLBOX IN PATH

try
    dtst = dataset(rand(10,100));       %attempt to create a dataset object
    props = properties(dtst);           %request the properties of said object
catch err
    if strcmp(err.identifier,'MATLAB:UndefinedFunction') %if dataset function not available then neither PLS Toolbox or Stats Toolbox are installed.
        

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
            sprintf('Cannot appropriately rejig path. Please manually place PLSToolbox above Stats Toolbox in path.');
            %fclose(fid);
            path(original_path);
            return
        end
        
    else    % If no stats toolbox entries have been found in path there is a more serious problem.
        sprintf('PLS Toolbox not on path. Please Install and try again.');
        %fclose(fid);
        path(original_path);
        return
    end
    
else
    sprintf('PLS Toolbox dataset objects are available. Continuing.')
    original_path = path;
end



%% LOAD PLS DATASET OBJECT

[dtst, ~, ~] = autoimport(dso_xml, 'xml');


%% EXTRACT THE MZ AXIS
ax = dtst.axisscale{dimension,1};%this assumes that the desired axis is the first axisscale set in the given dimension, it is possible to have multiple axisscales per dimension!

if ~isempty(ax)
    ax = ax(dtst.include{dimension,1});
    if isempty(ax)
        sprintf('no samples or variables included!')
        return
    end 
else
    sprintf('no axisscale in that dimension!')
    return
end


 

%% SAVE axis TO CSV
if dimension==1
    ax = ax'; %set to column if sample axis is desired
end

dlmwrite(outfile, ax, 'delimiter', ',', 'precision', 9)

%% RETURN PATH TO NORMAL
path(original_path);
rehash pathreset;
rehash toolbox;

return
end


