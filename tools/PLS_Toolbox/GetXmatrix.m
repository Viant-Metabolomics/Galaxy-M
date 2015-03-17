function GetXmatrix(dso_xml,outfile)
% -------------------------------------------------------------------
% GetXmatrix(dso_xml,outfile)
% 
% This is a very simple function created for use with the Birmingham Metabolomics Workflow in Galaxy. 
%
% inputs:
%	dso_xml	:	the path to an xml format PLS TOolbox dataset object as created by autoexport
%	outfile	:	the path to a csv file for output of the data matrix
% -------------------------------------------------------------------

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


%% EXTRACT THE X MATRIX
x = dtst.data(dtst.include{1,1},dtst.include{2,1});

%% SAVE X TO CSV
dlmwrite(outfile, x, 'delimiter', ',', 'precision', 9);

%% RETURN PATH TO NORMAL
path(original_path);
rehash pathreset;
rehash toolbox;

return
end


