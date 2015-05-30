function AddClasslist(dso_xml, listname, class_file, outfile)
% -------------------------------------------------------------------
% add_classlist_BB(dso_xml,listname,class_file,outfile, messagefile)
% 
% 
% -------------------------------------------------------------------

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

current_lists = dtst.classname(1,:);
n= length(current_lists);

%check that they aren't using a name that already exists (this could be better served with an 'overwrite' or 'reject' flag as input
for i = 1:n
	if strcmp(current_lists{i}, listname)
		disp('There already exists a class list with the name %s. Refusing to overwrite. Please choose a different name.',listname);
		path(original_path);
		return
	end
end

%add new classname
dtst.classname{1,n+1} = listname;

%% LOAD CSV FILE OF NEW CLASS LABELS
try
	classfid = fopen(class_file,'r');
	clist = textscan(classfid,'%s','Delimiter', ', ','MultipleDelimsAsOne',1);
	clist = clist{1}; %there ought to be one column of class labels... ought to be!
	fclose(classfid);
catch
	disp('Problem reading csv file');
	return
end


%% CHECK NUMBER OF CLASSES AND NUMBER OF ENTRIES
unique_classes = unique(clist);
sprintf('There are %i classes in this classlist', length(unique_classes));

if length(clist)~=size(dtst.data,1)
	sprintf('The number of class labels (%i) does not match the number of samples (%i). Quitting.',length(clist),size(dtst.data,1));
	return
end


%%CREATE THE CLASSLOOKUP ENTRY
class_lookup = cell(length(unique_classes),2);

for i = 1:length(unique_classes)
	class_lookup(i,:) = {i-1,unique_classes{i}};
end

dtst.classlookup{1,n+1} = class_lookup;	

%% CREATE THE NUMERICAL CLASS ID VECTOR
class_num = zeros(1,length(clist));

for i = 1:length(clist)
	for j = 1:length(unique_classes)
		if strcmp(clist{i},unique_classes(j))
			class_num(i) = j-1;
		end
	end
end
	
dtst.class{1,n+1} = class_num;
	


%% SAVE OUTPUT DATASET
autoexport(dtst, outfile, 'xml');

%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end


