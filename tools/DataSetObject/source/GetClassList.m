function GetClassList(dso_xml, listname, asnum, outfile)
% -------------------------------------------------------------------
% GetClassList(dso_xml,listname,asnum, outfile)
% 
% This is a function designed to work within the Birmingham Metabolomics Pipeline in Galaxy.  
%
% It extracts a classlist from the 1st (sample) dimension. This is intended for interoperability purposes - allowing the data to be passed to other pipelines that do not use the PLS Toolbox Dataset Object. 
%
% It makes use of the 'include' function of the dataset object. Changing the 'include' parameter for the Dataset will alter the output of this tool. 
%
% The output is a .csv file with one column and no headers.
%
% The option is given to have the class labels returned either using their text label e.g. 'control', 'control', 'QC'... or their integer value e.g. 1,1,8,...
%
%
% inputs:
%	dso_xml	:	the path to an xml format PLS Toolbox Dataset Object as output by autoexport()
%	listname:	a string to identify which classlist to return.
%	asnum	:	a boolean indicating the classes should be returned as text (0) or integer (1)
%	outfile	:	the path to a .csv file for outputting the two columns
% -------------------------------------------------------------------

if isa(asnum, 'char')
	asnum = str2num(asnum);
end

%% CHECK FOR PLS TOOLBOX IN PATH
if ~isdeployed
	try
	    dtst = dataset(rand(10,100));       %attempt to create a dataset object
	    props = properties(dtst);           %request the properties of said object
	catch err
	    if strcmp(err.identifier,'MATLAB:UndefinedFunction') %if dataset function not available then neither PLS Toolbox or Stats Toolbox are installed.
		

		disp('Matlab does not recognise dataset function. Neither PLSToolbox or Stats toolbox installed. Please amend and try again.');
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
		    disp('Cannot appropriately rejig path. Please manually place PLSToolbox above Stats Toolbox in path.');
		    path(original_path);
		    return
		end
		
	    else    % If no stats toolbox entries have been found in path there is a more serious problem.
		disp('PLS Toolbox not on path. Please Install and try again.');
		path(original_path);
		return
	    end
	    
	else
	    disp('PLS Toolbox dataset objects are available. Continuing.')
	    original_path = path;
	end
end


%% LOAD PLS DATASET OBJECT

[dtst, name, source] = autoimport(dso_xml, 'xml');

%% FIND THE CLASS LIST
classnames = dtst.classname(1,:); %only looking at the 1st dimension, sample classes.

matches = find(strcmp(listname, classnames)); %produces the indices of any matching class lists

if isempty(matches)
	disp('No classlists by that name, please try again! Quitting')
	return
elseif length(matches)>1
	disp('Multiple classlists by that name, please clean dataset object! Quitting')
	return
end

% EXTRACT CLASSES
if asnum %boolean positive means return integers
    classes = dtst.class{1,matches}'; %use transpose' to return as column
    classes = classes(dtst.include{1,1}); %exclude and 'not included' samples

    if isempty(classes)
        disp('No samples included! Quitting')
	return
    else
        csvwrite(outfile,classes);
    end

else % boolean negative means return class labels in text format

    classes = dtst.classid{1,matches}'; %returns the string values for each class (transpsed to a column).
    classes = classes(dtst.include{1,1}); %only return 'include' samples
 
    if isempty(classes)
        disp('No classes included! Quitting')
        return
    else
        fid = fopen(outfile ,'wt');
        for i = 1:length(classes)
            fprintf(fid, '%s\n',classes{i});
        end
        fclose(fid);
    end
    

end

%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end


