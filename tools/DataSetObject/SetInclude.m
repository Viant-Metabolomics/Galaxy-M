function set_include(dso_xml, dimension, include_file, outfile)

% -------------------------------------------------------------------
% set_include_BB(dso_xml, dimension, include_file, outfile)
% -------------------------------------------------------------------
dso_xml, dimension, include_file, outfile

if isa(dimension, 'char')
	dimension = str2num(dimension);
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

[dtst, ~, ~] = autoimport(dso_xml, 'xml');

inc = dtst.include{dimension};


%% LOAD CSV FILE OF NEW INCLUDE VALUES
try
	M = csvread(include_file);
catch
	disp('Problem reading csv file');
	return
end


%% COMPARE NEW INCLUDE TO EXISTING

if length(M)~=size(dtst.data,1)
	sprintf('The intended include list does not match the dimensions of the current data. Should have %i entries - 1 for include, 0 for exclude.',size(dtst.data,1));
	return
else
	dtst.include{dimension} = unique(find(M));
end




%% CREATE DATASET FOR OUTPUT

dtst_out=dtst;

%% SAVE OUTPUT DATASET
autoexport(dtst_out, outfile, 'xml');


%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end


