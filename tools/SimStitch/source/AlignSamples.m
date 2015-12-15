function AlignSample(fileList_in, html_indir, blankPPM, samplePPM, combPPM, out_pm, out_pm_row, out_pm_col)
%
%
% AlignSample_BB(fileList_in, html_indir,repfilter_file, blankPPM,samplePPM, combPPM,outfile)
%
%   A function that aligns the blanks with each other, the samples with
%   each other and then the blanks and the samples together. 
%
%   This is placed in the pipeline after ReplicateFilter and before
%   FlagBlankPeaks. 
%
%
%   Inputs - 
%       1, fileList_in: the xml from FileListManager
%       2, html_indir:	the full path to a directory, holding files referenced in html_infile - as used by Galaxy.
%       3, blankPPM: 	the ppm error to be used for the blanks (default = 2)
%       4, samplePPM: 	the ppm error to be used for the samples (default = 2)
%       5, combinedPPM: the ppm error to be used for aligning samples to
%                       blanks (default = 2)
%       6, out_pm:  	full path for peakmatrix .txt file.
%       7, out_pm_row:	full path to file for row headers for peak matrix
%       8, out_pm_col:	full path to file for column headers for peak matrix
%       9,out_params:	full path to file for summary of parameters output as text file. 
%
%       BLANKSTR: This is only changeable within the script. Blanks must be
%       labelled 'blank'in FileListManager!!!
%
%   Outputs: Various files are created but no matlab variables are returned by this function.
%
%   
%    RDavidson 31/01/13 (adapted from FlagBlankPeaks and SampleFilter)
%
%   Version 3.0
%   

%% %%%%%%%%%%Parameters
if isa(blankPPM, 'char')
	PARAMS.BLANK_MAXSPREAD_PPM = str2num(blankPPM); % the maximum spread (ppm) for aligning blanks
else
	PARAMS.BLANK_MAXSPREAD_PPM = blankPPM;
end

if isa(samplePPM, 'char')
	PARAMS.SAMPLE_MAXSPREAD_PPM = str2num(samplePPM); % the maximum spread (ppm) for aligning samples
else
	PARAMS.SAMPLE_MAXSPREAD_PPM = samplePPM;
end

if isa(combPPM, 'char')
	PARAMS.COMBINED_MAXSPREAD_PPM = str2num(combPPM); % the maximum spread (ppm) for aligning blanks with samples.
else
	PARAMS.COMBINED_MAXSPREAD_PPM = combPPM;
end

RFPLDir = [html_indir, filesep];

%one more parameter!
BLANKSTR = 'blank';
%% SORT PATH (MOVE STATS TOOLBOX ABOVE PLSTOOLBOX IF NECESSARY)
if ~isdeployed
	original_path = path; %save original path

	% Attempt a basic cluster, using stats toolbox. 
	try
		Y = pdist(rand(1,100)','cityblock');
		Z = linkage(Y,'complete');
		test = cluster(Z,'criterion','distance','cutoff',PARAMS.SAMPLE_MAXSPREAD_PPM);

	catch err   %check errors if fail
		if strcmp(err.identifier,'MATLAB:TooManyInputs')       %likely error 1
		    sprintf('Problem with too many inputs to cluster commands. Probably path order. Attempting temporary rejig')
		    
		    rem = original_path;
		    stats_path = '';
		    rem_path = '';
		    while true
		        [str,rem] = strtok(rem,pathsep);
		        if isempty(str)
		            break
		        elseif strfind(str,['stats',filesep,'stats']) %covers both '...\toolbox\stats\stats' and '...\toolbox\stats\statsdemos'
		            stats_path = [stats_path,str,pathsep];
		        elseif strfind(str,['stats',filesep,'classreg']) %covers the 3rd and final stats toolbox folder (as of 07/08/2013)
		            stats_path = [stats_path,str,pathsep];
		        else
		            rem_path = [rem_path,str,pathsep];
		        end
		    end
		    
		    if ~isempty(stats_path) %check for no stats toolbox installed
		        path(stats_path,rem_path);
		        rehash pathreset;
		        rehash toolbox;
		        
		        try
		            Y = pdist(rand(1,100)','cityblock');
		            Z = linkage(Y,'complete');
		            test = cluster(Z,'criterion','distance','cutoff',PARAMS.SAMPLE_MAXSPREAD_PPM);
		            sprintf('success!')
		        catch err
		            sprintf('cannot fix problems with Cluster function. Possible path not resetting. Fail.')
		            return
		        end
		    else
		        sprintf('stats toolbox apparently not installed! Fail.')
		        return
		    end
		        
		elseif strcmp(err.identifier,'MATLAB:UndefinedFunction') % likely error 2
		    sprintf('Matlab does not recognise cluster functions. Either PLSToolbox or Stats toolbox installed. Please amend and try again.')
		    return
		end
		
	end
end


%% %%%%%%%%%%Load Data (from replicate filter output)
fileList = ImportFileListXML(fileList_in);

fileList.filtPkDir = [html_indir,filesep];

%load spectra from replicate filter output folder
count = 0;

for si=1:fileList.nReplicates:length(fileList.Samples)
    count = count+1;
    % create filename 
    fname = [];
    
    for fi=si:(si+fileList.nReplicates-1)
        try
            fname = [fname,'_',fileList.Samples(1,fi).ID];
        catch
            fprintf('Error constructing filename: is batchList.nReplicates set correctly?\n');
            rethrow(lasterror);
        end
    end
    fname = ['RFPL_',fname,'.txt'];
    
    %open file in replicate filter output folder
    fid = fopen([RFPLDir,fname],'r');
    temp = textscan(fid,'%f\t%f\t%d\r\n','headerLines',1);
    fclose(fid);
    if si==29
        si;
    end
    %store spectral data for alignment
    sampleGpPks.name{count} = fname;
    sampleGpPks.blankflag{count} = strncmpi(fileList.Samples(1,si).ID,BLANKSTR,length(BLANKSTR));
    sampleGpPks.mz{count} = temp{1}.';
    sampleGpPks.y{count} = temp{2}.';
    sampleGpPks.n{count} = temp{3}.';
    
end

blank_index = find([sampleGpPks.blankflag{:}]);
sample_index = find([sampleGpPks.blankflag{:}]==0);

Num_blanks = 1;
Num_samples = 1; %assuming that there are both blanks and samples.

%% %%%%%%%%%%%%%%%%%%%%%%%%%align all blanks
if length(blank_index)>1
    Ci_blanks = PeakListCluster_0_3(sampleGpPks.mz(blank_index), PARAMS.BLANK_MAXSPREAD_PPM);
    
    
    fOptions_blanks.mFilt = false; % remove peaks occurring multiple times in same spectrum
    fOptions_blanks.nFilt = false; % peak must be present in n spectra
    fOptions_blanks.remCl = false;  % remove clusters filtered (no filtering anyway)
    fOptions_blanks.sort = true;   % sort output by MZbar
    [MZ_blanks, MZbar_blanks, Y_blanks, Ybar_blanks, N_blanks] = Clusters2Matrices_0_2(Ci_blanks, sampleGpPks.mz(blank_index), sampleGpPks.y(blank_index), [], fOptions_blanks);
    for i = 1:size(Y_blanks,1)
        Y_blanks_cell{i}= Y_blanks(i,:);
    end
elseif length(blank_index)==1
    MZ_blanks = sampleGpPks.mz{blank_index};
    MZbar_blanks = MZ_blanks;
    Y_blanks = sampleGpPks.y{blank_index};
    Ybar_blanks = Y_blanks;
    Y_blanks_cell{1}=Y_blanks;
else
    Num_blanks = 0;%no blanks...
end
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%align all samples
if length(sample_index)>1
Ci_samples = PeakListCluster_0_3(sampleGpPks.mz(sample_index), PARAMS.SAMPLE_MAXSPREAD_PPM);


fOptions_samples.mFilt = false; % remove peaks occurring multiple times in same spectrum
fOptions_samples.nFilt = false; % peak must be present in n spectra
fOptions_samples.remCl = false;  % remove clusters filtered (no filtering anyway)
fOptions_samples.sort = true;   % sort output by MZbar
[MZ_samples, MZbar_samples, Y_samples, Ybar_samples, N_samples] = Clusters2Matrices_0_2(Ci_samples, sampleGpPks.mz(sample_index), sampleGpPks.y(sample_index), [], fOptions_samples);

for i = 1:size(Y_samples,1)
    Y_samples_cell{i}= Y_samples(i,:);
end

elseif length(sample_index)==1
   MZ_samples = sampleGpPks.mz{sample_index};
   MZbar_samples = MZ_samples;
   Y_samples = sampleGpPks.y{sample_index};
   Ybar_samples = Y_samples;
   Y_samples_cell{1} = Y_samples;
else
    Num_samples=0; %no samples!
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%align samples and blanks

if Num_blanks & Num_samples %ie if there are both types of spectrum
mz_mixed = {MZbar_blanks, MZbar_samples};
y_mixed = {Ybar_blanks, Ybar_samples};
Ci_mixed = PeakListCluster_0_3(mz_mixed, PARAMS.COMBINED_MAXSPREAD_PPM);

%create a cluster index (ci), mz list and intensity list for each sample,
%but in their original order
Ci_mixed_expand = cell(1,size(sampleGpPks.name,2));
mz_mixed_expand = Ci_mixed_expand;
y_mixed_expand = Ci_mixed_expand;

Ci_mixed_expand(blank_index) = repmat(Ci_mixed(1),1,length(blank_index));
Ci_mixed_expand(sample_index) = repmat(Ci_mixed(2),1,length(sample_index));
mz_mixed_expand(blank_index) = repmat(mz_mixed(1),1,length(blank_index));
mz_mixed_expand(sample_index) = repmat(mz_mixed(2),1,length(sample_index));
y_mixed_expand(blank_index) = Y_blanks_cell;
y_mixed_expand(sample_index) = Y_samples_cell;

%combine the aligned blanks and samples! 
fOptions_mixed.mFilt = true;  % remove overlapping peaks
fOptions_mixed.nFilt = false; % don't filter by peak in minimum number of groups
fOptions_mixed.remCl = false; % don't remove clusters - set N to zero if overlapping
fOptions_mixed.sort = true;  % don't sort so output corresponds to input still
[MZ_mixed_final, MZbar_mixed_final, Y_mixed_final, Ybar_mixed_final, N_mixed_final] = Clusters2Matrices_0_2(Ci_mixed_expand, mz_mixed_expand, y_mixed_expand, [], fOptions_mixed);

else
    %if either blanks or samples are missing the final MZ and Y values have
    %already been calculated.
    if Num_samples
        MZ_mixed_final = MZ_samples;
        MZbar_mixed_final = MZbar_samples;
        Y_mixed_final = Y_samples;
        Ybar_mixed_final = Ybar_samples;
    else %must mean there are blanks (but no samples)(which is a little odd)
        MZ_mixed_final = MZ_blanks;
        MZbar_mixed_final = MZbar_blanks;
        Y_mixed_final = Y_blanks;
        Ybar_mixed_final = Ybar_blanks;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save the output to sample_filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%output folder


% Generate files
fprintf('Saving output\n');



%column headers
fid = fopen(out_pm_col,'w');
if ~fid, error('Cannot create file'); end
if iscell(MZbar_mixed_final) %this happens if there is only one sample or blank
    MZbar_mixed_final = cell2mat(MZbar_mixed_final);
end
fprintf(fid,'%0.5f\r\n',MZbar_mixed_final);
fclose(fid);
fprintf('done\n');

%table
fid = fopen(out_pm,'wt');
if ~fid, error('Cannot create file'); end
if iscell(Y_mixed_final)
    Y_mixed_final = cell2mat(Y_mixed_final);
end
for si=1:size(Y_mixed_final,1)
    fc = fprintf(fid,'%0.7g\t',Y_mixed_final(si,:));
    fc = fc & fprintf(fid,'\n');
    if ~fc, error('Cannot write to file'); end
end

fclose(fid);
fprintf('done\n');

%row headers
fid = fopen(out_pm_row,'wt');
if ~fid, error('Cannot create file'); end
blank_cell = {'sample', 'blank'};
for si=1:length(sampleGpPks.name)
    fc = fprintf(fid,'%s\t',blank_cell{sampleGpPks.blankflag{si}+1});
    fc = fc & fprintf(fid,'%s\n',sampleGpPks.name{si});
    if ~fc, error('Cannot write to file'); end
end

fclose(fid);
fprintf('done\n');

%parameters
%fid = fopen(out_params,'w');
%if ~fid, error('Cannot create file'); end
%fprintf(fid,'FILE VERSION:\t%s',mfilename); 
%fprintf(fid,'\r\n\r\nPARAMS:\r\n');
%fn = fieldnames(PARAMS);
%for i=1:length(fn)
%    val = PARAMS.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
%fclose(fid);
%fprintf('done\n');

% % Update fileListmanager paths
% xDoc = xmlread(fileList_in);
% xDoc.getElementsByTagName('ALPMDir').item(0).getChildNodes.setTextContent(batchList.ASDir);
% write_XML_no_extra_lines(fileList_in,xDoc);


%% RETURN PATH TO ORIGINAL STATE (In case of rejigging due to cluster function issues)
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end


%% FINISHED

fprintf('Finished.\n');
return
end

