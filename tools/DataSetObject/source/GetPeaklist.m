function GetPeaklist(dso_xml, outfile)
% -------------------------------------------------------------------
% GetPeaklist(dso_xml,outfile)
% 
% This is a function designed to work within the Birmingham Metabolomics Pipeline in Galaxy. 
%
% This output is intended to be used as an input for the MI Pack peak identification tools. 
%
% It extracts the axisscale in the 2nd dimension, intended to be the MZ values of your mass spectra, and also extracts the average intensity of each peak across all 'included' data. 
%
% Thus it is possible to manipulate the number of entries in the output by changing which peak variables are 'included', and the average intensities by altering which samples are 'included', where 'include' is a function applicable to PLS Dataset Objects. 
%
% The output is a tab delimited file with two columns, one for MZ value and one for Intensity. These text headers are included in the file.
%
%
% inputs:
%	dso_xml	:	the path to an xml format PLS Toolbox Dataset Object as output by autoexport()
%	outfile	:	the path to a .tab file for outputting the two columns
%
% -------------------------------------------------------------------

%fid = fopen(messagefile,'a');
%fprintf(fid, 'Messages for get_peaklist:');

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
'auto'
end


%% LOAD PLS DATASET OBJECT
[dtst, name, source] = autoimport(dso_xml, 'xml');

%% EXTRACT THE MZ AXIS
mz = dtst.axisscale{2,1}(dtst.include{2,1}); %this assumes that the MZ axis is the first axisscale set in the 2nd dimension, it is possible to have multiple axisscales per dimension!
intensity = mean(dtst.data(dtst.include{1,1},dtst.include{2,1}),1);
%% SAVE MZ TO CSV
fid_mz = fopen(outfile,'a');

fprintf(fid_mz, 'MZ\t Intensity\n'); %This header is necessary for interaction with MIPACK (Birmingham MEtabolomics Pipeline)
fprintf(fid_mz, '%f\t %f\n',[mz',intensity']');
fclose(fid_mz);


%Close the message output file
%fclose(fid);
 

%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end


