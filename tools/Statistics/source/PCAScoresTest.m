function PCAScoresTest(dso_xml, classlist, outfile_l, outfile_s, outfile_v, results, ncomp)
%
%  PCAScoresTestFree(dso_xml, classlist, outfile_l, outfile_s, outfile_v, results, ncomp)
%
%  inputs:
%	dso_xml:	Full path to an xml file containing a PLS Toolbox dataset object
%	classlist:	the name of the list of classes to use for PCA scores testing
%	outfile_l:	full path to a csv file to be used for storing loadings vectors
%	outfile_s:  full path to a csv file to be used for storing PC scores
%	outfile_v:  full path to a csv file to be used for storing Eigenvalues
%	results:	Full path to a .txt file for writing the scores test results
%	ncomp:		The number of components to use in the model. 
%
%  outputs:
%       The script does not return any Matlab variables but does write to the files specified by outfile, results and messages. 
%
%  PCA is performed with mean centring. Scaling is not currently applied here. It is assumed that GLOG transformation will have been applied beforehand.  
%  For 2 classes of samples, a t-test is performed for each PC giving a p value that describes the separation of the two classes for that PC.
%  For 3 or more classes of samples, ANOVA results are given for each PC and Tukey post-hoc testing is also reported.
%
%
%  J.Byrne, J.Kirwan A.Southam, M.Viant. 
%  School of Biosciences, University of Birmingham, August 2010
%
%
%  Adapted for Galaxy Project by R.L.Davidson 02/09/2013

if isa(ncomp, 'char')
	ncomp = str2num(ncomp);
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



% open output files for status reporting
%fid = fopen(messages,'a');
sprintf( 'Messages pertaining to running of PCA - Automatic with scores test');
fid_results = fopen(results, 'a');
fprintf(fid_results, 'Results pertaining to PCA- Automatic with scores test.\n');

%% CREATE PCA MODEL

%only include desired samples (inclusion is stored as a variable in PLS
%dataset object)
samples_included=dso.include{1,1}; %Again, for simplicity we will assume 1 include array for samples and one for peaks/variables so we can hard code {1,1}. 
x_mat = dso.data(samples_included,:);

%x_mat = x_mat + rand(size(x_mat)) % test
x_mat_MC = x_mat - repmat(mean(x_mat),size(x_mat,1),1);

[U,D,coeff] = svd(x_mat_MC,'econ');
score = U*D;
latent = diag(D'*D)/(size(x_mat_MC,1)-1);

% [coeff,score,latent] = princomp(x_mat);
% reduce based on desired number of PCs

cum_thresh = 0.7;
single_thresh = 0.15;

if ncomp==0 %user has asked script to automatically choose number of components
    cum_varz = cumsum(latent)./sum(latent);
    single_props = latent./(sum(latent));
    
    potential_pcs = intersect(find(cum_varz>cum_thresh),find(single_props>single_thresh)); %this is combining two rules of thumb for picking components
    if length(potential_pcs)==0
        lvs = min(find(cum_varz>cum_thresh)); %if there are no PCs after the cumulative variance threshold that have high individual variance, just use the cumulative variance threshold
    else
        lvs = max(potential_pcs); % if there are PCs with relatively high variance, after cumulative variance has reached a threshold, keep them in.
    end
    
elseif ncomp>=(size(dso.data,1)-1) %user has chosen too many PCs 
	lvs = size(dso.data,1)-1; %set number of components to maximum
	sprintf('Too many PCs requested by user. Choosing full model with %i PCs.\n',lvs);
else
	lvs= ncomp; %user has specified a legitimate number of PCs so we shall use that value.
end


% reduce outputs to match number of components
loadings = coeff(:,1:lvs)';
scores = score(:,1:lvs);
values = latent(1:lvs);

% add MZ values to loadings output
ax = dso.axisscale{2,1};
ax = ax(dso.include{2,1});

loadings = [ax;loadings];

%% SAVE OUTPUT MODEL

%save the 3 outputs as 3 separate CSV files. 

csvwrite(outfile_l,loadings);
csvwrite(outfile_s,scores);
csvwrite(outfile_v,values);




%% PERFORM PCA SCORES TEST

fprintf(fid_results,'%s \n \n', '**************   PCA SCORES TEST RESULTS   **************');

nsamps=size(dso.data,1);

%Select classlist name as specified by user

%all_classlists = model.detail.classname(1,:); %2nd dimension is used for variable classes which are not an option here.
all_classlists = dso.classname(1,:);
classlist_index = 0;
for i = 1:length(all_classlists)
    if strcmp(all_classlists{i},classlist)
        classlist_index=i;% This will select the last instance of the classlist with that name.  
         
    end
end

if classlist_index %in the case where an index was found
    
    classes=dso.class{1,classlist_index}; %again, 1st dimension is for sample classes and is only one being used here.
    class_lookup=dso.classlookup{1,classlist_index};
    labels=dso.label{1,1};  %These are sample names, there could be multiple labels per sample but for simplicity we're only going to use the first set:  label{1,1} 
    
    if ((length(classes))==0)    %  This script uses the class information from the input model
        sprintf('ERROR:  No class information defined in the classlist. %s\n', classlist);% appears in red
    end
    

    classes=classes(samples_included);
    [sorted_classes,idx]=sort(classes);
    sorted_scores=scores(idx,:);
    [n_included_samps,n_pcs]=size(sorted_scores);
    unique_class_ids=unique(sorted_classes);
    no_unique_class_ids=length(unique_class_ids);   %The number of unique class ID's.
    n_excluded=nsamps-n_included_samps;
    if (n_excluded==0)
    elseif (n_excluded==1)
        fprintf(fid_results,'%s %s %s \n', 'The data set consists of', num2str(nsamps),'samples of which 1 was excluded from the model.');
    else
        fprintf(fid_results,'%s %s %s %s %s \n', 'The data set consists of', num2str(nsamps),'samples of which', num2str(n_excluded), 'were excluded from the model.');
    end
    
    if (n_excluded==0)
        fprintf(fid_results,'%s \n', 'All samples were included in the model.');
    else  %  we have some exclusions of samples
        
        all_samples_vector=1:nsamps;
        all_samples_vector(samples_included)=0;  % this leaves only the excluded samples having non-zero entries.
        excluded_samples=find(all_samples_vector);
        
        if (strcmp(dso.label{1,1},''))  %  there are no sample names.  Just report the initial positions of the excluded samples.
            fprintf(fid_results,'%s \n', 'No sample names were provided.');
            fprintf(fid_results,'%s ', 'Initial row positions of excluded samples:');
            fprintf(fid_results,'%s', num2str(excluded_samples));
            fprintf(fid_results,'\n');
            
        else   %  There are sample names available.
            
            fprintf(fid_results,'%s', 'Excluded samples: ');
            
            for i=1:n_excluded
                fprintf(fid_results,'%s',labels(excluded_samples(i),:),'; ');
            end
            fprintf(fid_results,'\n');
        end
    end
    fprintf(fid_results,'%s %s %s\n', 'The model was constructed using data from',num2str(n_included_samps),'samples.');
    fprintf(fid_results,'%s %s %s', 'The samples were grouped into',num2str(no_unique_class_ids),'classes: ');
    for i=unique_class_ids
        fprintf(fid_results, '%s ',class_lookup{find([class_lookup{:,1}]==i),2},';' );
    end
    fprintf(fid_results,'\n');
    
    %   GET THE SAMPLE NAMES from dso.label
    if (strcmp(dso.label{1,1},''))   % i.e. there are no sample names stored
    else  % WE HAVE SAMPLE NAMES
        %  reorder the samples for each unique class, find the samples of that class.
        labels=labels(samples_included,:);
        sorted_labels=labels(idx,:);
        
        fprintf(fid_results,'\n%s \n', 'CLASS MEMBERSHIP DETAILS:');
        for i=unique_class_ids
            members_of_class=find(sorted_classes==i);
            fprintf(fid_results, '%s %s ', class_lookup{find([class_lookup{:,1}]==i),2},'(',num2str(length(members_of_class)),'of)',':');
            
            for j=1:length(members_of_class)
                fprintf(fid_results, '%s ',sorted_labels(members_of_class(j),:),';');
            end
            fprintf(fid_results, '\n');
        end
        
    end
    
    %  DISPLAY TOTAL VARIANCE
    cum_var = sum(latent(1:lvs))/sum(latent);

    fprintf(fid_results,'\n%s %0.2f%% %s %s %s\n', 'The tested model explains',cum_var, 'of the total variance using', num2str(n_pcs),'PCs.');
    %
    %    CHOOSE WHICH TYPE OF TESTING TO PERFORM BASED ON THE NUMBER OF CLASSES
    %    INCLUDED IN THE MODEL.  IF 1, NO TESTING, IF 2, TTEST, IF 3 ANOVA AND
    %    POST-HOC TUKEY KRAMER.
    
    if (no_unique_class_ids==1)
        sprintf('There is only one sample class in this model.  Is that sensible?  No statistics performed.\n')% appears in red
        %results_structure=[];
        
    elseif (no_unique_class_ids==2)    %DO TTEST2
        
        fprintf(fid_results,'\n%s %s %s \n', 'T-test p-value results for the',num2str(n_pcs), 'PC PCA model:');
        
        for i=1:2
            nidx=length(find(sorted_classes==unique_class_ids(i)));
            if(i==1)
                class1_sample_ids(1:nidx,1)=find(sorted_classes==unique_class_ids(i));
            elseif(i==2)
                class2_sample_ids(1:nidx,1)=find(sorted_classes==unique_class_ids(i));
            end
        end
        
        for i=1:n_pcs      %For each PC
            class_1_scores_vector=sorted_scores(class1_sample_ids,i);
            class_2_scores_vector=sorted_scores(class2_sample_ids,i);
            [~, p]=ttest2(class_1_scores_vector,class_2_scores_vector);
            %results_structure.ttest.pcs(i).h=h;
            %results_structure.ttest.pcs(i).p=p;
            var_captured=latent(i)/sum(latent);
            fprintf(fid_results,'%s %s %s %s %s %0.2f%% %s\n', 'PC', num2str(i),' p-value = ',num2str(p),'(', var_captured,'of variance)');
        end
        
    elseif (no_unique_class_ids>=3)  %%%   DO AN ANOVA ON THE DATA
        
        %results_structure.anova.results_overview={};
        %results_structure.anova=[];
        
        %clear options
        %options.display = 'off';
        %options.plots = 'none';
        %options.preprocessing = {'mean center' 'mean center'};
        %options.algorithm = 'sim';
        %options.discrim = 'no';
        %options.structureoutput = 'yes';
        
        %%%  Create the ID matrix, which contains columns of element ID's which are not necessarily the same number of rows long.  These ID's
        %%%  are the addresses in the sorted_scores matrix where the scores from each of the different classes is found
        
        fprintf(fid_results,'\n%s %s %s\n', 'ANOVA performed.  Results for the',num2str(n_pcs), 'PC PCA model:');
        
        for i=1:n_pcs      %For each PC
            
            [pvalue,table,stats] = anova1(sorted_scores(:,i),sorted_classes,'off');
            %results_structure.anova.pcs(i).anova1_pvalue=pvalue;
            %results_structure.anova.pcs(i).anova1_table=table;
            %results_structure.anova.pcs(i).anova1_stats=stats;
            [comparison,~] = multcompare(stats, 'ctype', 'tukey-kramer','display','off');
            
            comparison_output=comparison;
            comparison_output(:,1)=unique_class_ids(comparison_output(:,1))';
            comparison_output(:,2)=unique_class_ids(comparison_output(:,2))';
            
            [nrows ~]=size(comparison_output);
            significants={};
            
            
            for j=1:nrows
                significants{j}='';
                if((comparison_output(j,3)*comparison_output(j,5))>0)  %  same signs will give a +ve result, different signs will be -ve.
                    significants{j}='significant';
                end
            end
            
            %results_structure.anova.pcs(i).multcompare_comparison=comparison;
            %esults_structure.anova.pcs(i).multcompare_means=means;
            
            var_captured=latent(i)/sum(latent);
            fprintf(fid_results,'%s\n', '_________________________________________________________');
            fprintf(fid_results,'\n%s %s %s %0.2f%% %s\n', 'Results for PC',num2str(i),'(',var_captured,'of variance)');
            fprintf(fid_results,'%s %s\n', 'p-value = ',num2str(pvalue));
            f=sprintf('%s',table{2,5});
            fprintf(fid_results,'%s %s \n', 'F-value = ', f);
            
            fprintf(fid_results,'%s \n \n', 'Tukey-Kramer post-hoc analysis performed.  Comparison table:');
            fprintf(fid_results,[repmat('%8s \t',1,5),'\n'],'group ',' group ','lower ','estimated ','upper ');
            
            for j=1:nrows
                class1 = class_lookup{find([class_lookup{:,1}]==comparison_output(j,1)),2};
                class2 = class_lookup{find([class_lookup{:,1}]==comparison_output(j,2)),2};
                
                fprintf(fid_results,['%8s \t %8s \t' ,repmat('%8.4f \t',1,3),'%s \n'],class1,class2,comparison_output(j,3:5),significants{j});
            end
            
        end
    end

    
else % where the classlist name does not match any classlist label
    sprintf('ERROR:  No class information defined in the classlist. %s\n', classlist);% appears in red
    fprintf(fid_results, 'No scores test results due to mismatch with classlist name');
end




%make an output for Galaxy
fclose(fid_results);
%fclose(fid);

%% RETURN PATH TO NORMAL
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end
return

end
