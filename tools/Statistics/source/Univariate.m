function [h, p, crit_p, adj_p, fold_change] = Univariate(dso_xml, fdr, controlclass, filename, parametric)
%
%   This function applies either a t-test (two-tailed) or anova (1-way)
%   depending on the number of classes present.%
%
%   The Benjamini-Hochberg correction is applied at a false discovery rate
%   set by the user. It is not necessary to set a confidence value for the
%   t-test or ANOVA.
%
%   Fold changes are reported with reference to the control class. Therefore
%   there will be 1 row less than the total number of classes.
%
%       NB: the t-test performed is an unpaired, 2-tailed t-test. It
%       assumes that the variance of the two classes is equal.
%
%       NB: The anova is a 1-way ANOVA that is best when the classes are of
%       equal size. 
%
%   inputs:
%       dso_xml = a PLS toolbox dataset object file (xml format)
%       fdr = the false discovery rate you wish to apply e.g. 0.05, 0.1 etc
%       controlclass = a string containing the name of the class to be used
%                       as control group (for fold change).
%                       e.g. controlclass = 'Class 1';
%       filename = a string containing the name of the output file
%       parametric =    set to 1 where data is normal (t-test or ANOVA)
%                       set to 0 where data is non-normal (ranksum or
%                       kruskalwallis)
%                       DEFAULT = 1;
%
%   outputs:
%       h = a vector of 1s and 0s indicating which peaks have passed (1) or
%           failed (0) the Benjamini-Hochberg correction
%       p = the original, unadjusted P values for each peak
%       crit_p = the adjusted P threshold. 
%               Therefore, significant peaks can be calculated by using
%               sig_peak_index = find(p<=crit_p);
%       adj_p = P values for each peak, adjusted so that 0.05 is a valid
%               cut-off. Therefore, significant peaks can be calulated
%               using:
%               sig_peak_index = find(adj_p<=0.05);
%       fold_change =   For t-tests, a single vector representing the fold
%                       change from control to test ( >1, increase in test; <1, decrease in test)
%                       For ANOVA, there will be columns for each test class,
%                       ordered by class label, top to bottom 
%                       e.g. classes = [2,2,2,3,3,3,1,1,1], control = 3,
%                       then fold change will have a row for class 1 and
%                       then a row for class 2. 
%       output file =   this is a tab delimited file that can be read by
%                       MS Excel. The top row will have headers describing
%                       each column. They are Peak ID, P, H, ADJ_P and a
%                       column for each fold-change. 
%

if isa(fdr, 'char')
	fdr = str2num(fdr);
end

if isa(parametric, 'char')
	parametric = str2num(parametric);
end


if nargin <5
    parametric = 1;
    disp('no parametric option detected so using default assumption of normality')
end



%% CHECK FOR PLS TOOLBOX IN PATH
if ~isdeployed
	try
	    dtst = dataset(rand(10,100));       %attempt to create a dataset object
	    props = properties(dtst);           %request the properties of said object
	catch err
	    if strcmp(err.identifier,'MATLAB:UndefinedFunction') %if dataset function not available then neither PLS Toolbox or Stats Toolbox are installed.
		
		%fid = fopen(messages,'w');
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
		    %fid = fopen(messages,'w');
		    sprintf('Cannot appropriately rejig path. Please manually place PLSToolbox above Stats Toolbox in path.');
		    %fclose(fid);
		    path(original_path);
		    return
		end
		
	    else    % If no stats toolbox entries have been found in path there is a more serious problem.
		%fid = fopen(messages,'w');
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

%% GET DATA FROM FILELIST OBJECT
%load data
[dso, name, source] = autoimport(dso_xml, 'xml');

if isa(dso, 'dataset')
    store_dso = dso;
    peaks = dso.data;
    
    %include = dso.include{1};
    
    % Automatically detect QC and blank samples and exclude from the
    % analysis
    [n,~] = size(dso);
    QC_samples = [];
    Blank_samples = [];
    for i = 1:n
        if strcmp(dso.classid{1}(i),'QC') || strcmp(dso.classid{1}(i),'qc') || strcmp(dso.classid{1}(i),'Qc') || strcmp(dso.classid{1}(i),'qC')
            QC_samples = [QC_samples i];
        elseif strcmp(dso.classid{1}(i),'blank') || strcmp(dso.classid{1}(i),'BLANK')
            Blank_samples = [Blank_samples i];
        end
    end
    include = 1:n; include(sort([QC_samples,Blank_samples],'ascend')) = [];
    
    include_p = dso.include{2};
    
    if isempty(include_p)
        sprintf('no peaks included in dataset object!')
        h = [];
        p = [];
        crit_p = [];
        adj_p = [];
        fold_change = [];
        return;
    end
    
    if isempty(include)
        sprintf('no samples "included" in dataset object!')
        h = [];
        p = [];
        crit_p = [];
        adj_p = [];
        fold_change = [];
        return;
    end
    
    peaks = peaks(include,include_p);
    peaks(find(peaks==0))=NaN;
    
    classlist = [];
    control = [];
    for i = 1:length(dso.classlookup(1,:))
        if ~isempty(control)
            break
        else
            for j = 1:size(dso.classlookup{1,i},1)
                if strcmp(dso.classlookup{1,i}{j,2},controlclass)
                    classlist = i;
                    control = dso.classlookup{1,i}{j,1};
                    
                end
            end
        end
    end
    
    if ~isempty(classlist)
        classes = store_dso.class{1,classlist};
        classes = classes(include);
        class_labels = unique(classes);
        if length(class_labels)<2
            sprintf('Not enough classes in class list')
            h = [];
            p = [];
            crit_p = [];
            adj_p = [];
            fold_change = [];
            return;
        end
    else
        sprintf('Cannot find class with ''controlname'' in dataset object!')
        h = [];
        p = [];
        crit_p = [];
        adj_p = [];
        fold_change = [];
        return;
    end
    

end
    


q = fdr;



numClasses = length(class_labels);


if numClasses==2
    if parametric
        disp('performing t-test')
    else
        disp('performing wilcoxon rank sum test')
    end
    control_g = peaks(find(classes==control),:);
    test = peaks(find(classes~=control),:);
    class_labels(find(class_labels==control))=[]; %to help with writing out to file at end
    
    for i = 1:size(peaks,2)
        if parametric
            
            [H(i),p(i)] = ttest2(control_g(:,i),test(:,i));
        else
            
            [p(i),H(i)] = ranksum(control_g(:,i),test(:,i));
        end                     
        fold_change(i) = nanmean(test(:,i))/nanmean(control_g(:,i));
    end
    fold_change = fold_change';
    method = 'pdep';
    report = 'yes';
    [h, crit_p, adj_p]=fdr_bh(p,q,method,report);
else
    if parametric
        sprintf('More than two classes, using 1 way ANOVA');
    else
        sprintf('More than two classes, using KruskalWallis test');
    end
    
    for i = 1:numClasses
        ind{i} = find(classes==class_labels(i));
    end
    control_ind = find(class_labels==control)
    control_g = ind{control_ind};
    ind(control_ind) = [];
    class_labels(control_ind) = [];
       
    
    for j = 1:size(peaks,2)
        if parametric
            
            p(j) = anova1(peaks(:,j), classes, 'off');
        else
            p(j) = kruskalwallis(peaks(:,j),classes, 'off');
        end
        for k=1:length(class_labels)
            fold_change(j,k) = nanmean(peaks(ind{k},j))/nanmean(peaks(control_g,j));
        end 
    end
    method = 'pdep';
    report = 'yes';
    [h, crit_p, adj_p]=fdr_bh(p,q,method,report);
    
    
  
end
    
label=dso.axisscale{2};
if isempty(label)
    label = 1:size(dso.data,2);
end
label = label(include_p)';
    
h=h';
adj_p=adj_p';
p=p';
crit_p = crit_p';
   
outMat = [label,h,p,adj_p,fold_change];

fid = fopen(filename,'w');
if fid<0, error('Cannot create output file - Is file open already?'); 
else
    classnames = dso.classlookup{1,classlist}(:,2);
    for i = 1:size(fold_change,2)
        class_header{i} = classnames{find(cell2mat(dso.classlookup{1,classlist}(:,1))==class_labels(i))};
    end
    
    fout_headers =  'Peak ID\tH\tP\tADJ_P';
    fout_data =     '%.5f\t%i\t%0.2f\t%0.6f';
    for i = 1:length(class_header)
        fout_headers = [fout_headers,'\t',class_header{i}];
        
        fout_data = [fout_data,'\t%.6f'];
    end
    fout_headers = [fout_headers,'\n'];
    fout_data = [fout_data,'\n'];
    
    fprintf(fid,fout_headers);
    fprintf(fid,fout_data,outMat');

    fclose(fid);
    fprintf('done\n');
    
end

        
        
        
