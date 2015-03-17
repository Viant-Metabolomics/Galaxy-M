function AlignSample_BB(fileList_in,html_infile, html_indir, blankPPM,samplePPM, combPPM,out_pm, out_pm_row, out_pm_col)
%
%
% AlignSample_BB(fileList_in,html_infile, html_indir,repfilter_file, blankPPM,samplePPM, combPPM,outfile)
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
%       2, html_infile:	the full path to an html file output by Rep Filter	
%       3, html_indir:	the full path to a directory, holding files referenced in html_infile - as used by Galaxy.
%       4, blankPPM: 	the ppm error to be used for the blanks (default = 2)
%       5, samplePPM: 	the ppm error to be used for the samples (default = 2)
%       6, combinedPPM: the ppm error to be used for aligning samples to
%                       blanks (default = 2)
%       7, out_pm:  	full path for peakmatrix .txt file.
%       8, out_pm_row:	full path to file for row headers for peak matrix
%       9, out_pm_col:	full path to file for column headers for peak matrix
%       10,out_params:	full path to file for summary of parameters output as text file. 
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

PARAMS.BLANK_MAXSPREAD_PPM = blankPPM; % the maximum spread (ppm) for aligning blanks
PARAMS.SAMPLE_MAXSPREAD_PPM = samplePPM; % the maximum spread (ppm) for aligning samples
PARAMS.COMBINED_MAXSPREAD_PPM = combPPM; % the maximum spread (ppm) for aligning blanks with samples.

%one more parameter!
BLANKSTR = 'blank';
%% SORT PATH (MOVE STATS TOOLBOX ABOVE PLSTOOLBOX IF NECESSARY)

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



%% %%%%%%%%%%Load Data (from replicate filter output)
fileList = import_filelist_xml(fileList_in);

fileList.filtPkDir = [html_indir,filesep];

PARAMS.NUMREPS = fileList.numReps;     %the number of replicates per sample

%load spectra from replicate filter output folder
count = 0;
for si=1:PARAMS.NUMREPS:length(fileList.rawShort)
    count = count+1;
    % create filename 
    fname = fileList.setUID;
    
    for fi=si:(si+PARAMS.NUMREPS-1)
        try
            fname = [fname,'_',fileList.rawShort{fi}];
        catch
            fprintf('Error constructing filename: is PARAMS.NUMREPS set correctly?\n');
            rethrow(lasterror);
        end
    end
    fname = ['filteredPeaks_',fname,'.txt'];
    
    %open file in replicate filter output folder
    fid = fopen([fileList.filtPkDir,fname],'r');
    temp = textscan(fid,'%f\t%f\t%d\r\n','headerLines',1);
    fclose(fid);
    if si==29
        si;
    end
    %store spectral data for alignment
    sampleGpPks.name{count} = fname;
    sampleGpPks.blankflag{count} = strncmpi(fileList.spec(si).sampleID,BLANKSTR,length(BLANKSTR));
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


%% RETURN PATH TO ORIGINAL STATE (In case of rejigging due to cluster function issues)

path(original_path);
rehash pathreset;
rehash toolbox;




%% FINISHED

fprintf('Finished.\n');
return
end


function fileList_struct = import_filelist_xml(file_location)
% 
% function fileList_struct = import_filelist_xml(file_location)
%
% function to allow SimStitch pipeline to read/import xml representations of individual fileList entries 
%
% specifically aimed at the Galaxy Project implementation of SimStitch etc.
% inputs:
%   file_location:  the location of a fileList.xml file as produced by
%                   FileListManagerGUI
%
% outputs:
%   fileList_struct:a matlab struct conforming to SimStitch's fileList
%                   structure
%
% R.L.Davidson 
% 29/05/2014



try
    tree = xmlread(file_location);
catch
    errordlg('Failed to read XML file.','Import error');
end

if ~tree.hasChildNodes()
    errordlg('XML tree empty!','Import error');
    return
end

% this struct declaration is not very concise, but helps me remember the structure (could be
% declared like fileListStruct.spec is e.g.
% struct('setIdentifier',{},'rootIdentifier',{},...)
fileListStruct = [];
fileListStruct.setIdentifier = '';
fileListStruct.rootDirectory= '';
fileListStruct.rawDirectory= [];
fileListStruct.spec= struct('ID',{},'rawFile',{},'sampleID',{});
fileListStruct.analysisSplit= [];
fileListStruct.datDirectory= [];
fileListStruct.avgdTransDirectory= '';
fileListStruct.overlappingSpecDirectory= '';
fileListStruct.stitchedSpecDirectory= '';
fileListStruct.peakListsDirectory= '';
fileListStruct.filteredPeaksDirectory= '';
fileListStruct.blankFlaggedDirectory= '';
fileListStruct.combinedPeaksDirectory= '';
fileListStruct.internalCalFile= [];
fileListStruct.numReps= [];
fileListStruct.count= [];
fileListStruct.notes= [];
fileListStruct.spectraCount= [];
fileListStruct.DSO_Directory= '';
fileListStruct.MultivarModel_Directory= '';
fileListStruct.MetID_Directory= '';
fileListStruct.specHeaders= {};
fileListStruct.specData= {};


%i mislabeled the struct above, will rename for ease now (lazy!)
fileList_struct = fileListStruct;

%start to interrogate the xml input
if strcmp('fileList',tree.getChildNodes.item(0).getNodeName)
    fileList = tree.getChildNodes.item(0);
else
    errordlg('XML Document root is not fileList - stopping')
    return
end

try
    
    if ~isempty(fileList.getElementsByTagName('setIdentifier').item(0).getChildNodes.item(0))
	fileList_struct.setIdentifier = char(fileList.getElementsByTagName('setIdentifier').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('rootDirectory').item(0).getChildNodes.item(0))    
	fileList_struct.rootDirectory = char(fileList.getElementsByTagName('rootDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('rawDirectory').item(0).getChildNodes.item(0))
        fileList_struct.rawDirectory = char(fileList.getElementsByTagName('rawDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    fileList_struct.analysisSplit = str2num(fileList.getElementsByTagName('analysisSplit').item(0).getChildNodes.item(0).getData);
    
    if ~isempty(fileList.getElementsByTagName('datDirectory').item(0).getChildNodes.item(0))
        fileList_struct.datDirectory = char(fileList.getElementsByTagName('datDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('avgdTransDirectory').item(0).getChildNodes.item(0))
        fileList_struct.avgdTransDirectory = char(fileList.getElementsByTagName('avgdTransDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('overlappingSpecDirectory').item(0).getChildNodes.item(0))
        fileList_struct.overlappingSpecDirectory = char(fileList.getElementsByTagName('overlappingSpecDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('stitchedSpecDirectory').item(0).getChildNodes.item(0))
        fileList_struct.stitchedSpecDirectory = char(fileList.getElementsByTagName('stitchedSpecDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('peakListsDirectory').item(0).getChildNodes.item(0))
        fileList_struct.peakListsDirectory = char(fileList.getElementsByTagName('peakListsDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('filteredPeaksDirectory').item(0).getChildNodes.item(0))
        fileList_struct.filteredPeaksDirectory = char(fileList.getElementsByTagName('filteredPeaksDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('blankFlaggedDirectory').item(0).getChildNodes.item(0))
        fileList_struct.blankFlaggedDirectory = char(fileList.getElementsByTagName('blankFlaggedDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('combinedPeaksDirectory').item(0).getChildNodes.item(0))
        fileList_struct.combinedPeaksDirectory = char(fileList.getElementsByTagName('combinedPeaksDirectory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('internalCalFile').item(0).getChildNodes.item(0))
        fileList_struct.internalCalFile = char(fileList.getElementsByTagName('internalCalFile').item(0).getChildNodes.item(0).getData);
    end
    
    
    fileList_struct.numReps = str2num(fileList.getElementsByTagName('numReps').item(0).getChildNodes.item(0).getData);
    
    fileList_struct.count = str2num(fileList.getElementsByTagName('count').item(0).getChildNodes.item(0).getData);
    
    if ~isempty(fileList.getElementsByTagName('notes').item(0).getChildNodes.item(0))
        fileList_struct.notes = char(fileList.getElementsByTagName('notes').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('DSO_Directory').item(0).getChildNodes.item(0))
        fileList_struct.DSO_Directory = char(fileList.getElementsByTagName('DSO_Directory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('MetID_Directory').item(0).getChildNodes.item(0))
        fileList_struct.MetID_Directory = char(fileList.getElementsByTagName('MetID_Directory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('MultivarModel_Directory').item(0).getChildNodes.item(0))
        fileList_struct.MultivarModel_Directory = char(fileList.getElementsByTagName('MultivarModel_Directory').item(0).getChildNodes.item(0).getData);
    end
    
    if ~isempty(fileList.getElementsByTagName('specHeaders').item(0).getChildNodes.item(0))
        header_char = char(fileList.getElementsByTagName('specHeaders').item(0).getChildNodes.item(0).getData);
        specHeaders = {};
        while ~isempty(header_char)
            [specHeaders{end+1},header_char] = strtok(header_char,',');
        end
        fileList_struct.specHeaders = specHeaders;
    end
    
    if ~isempty(fileList.getElementsByTagName('specData').item(0).getChildNodes.item(0))
       data_char = char(fileList.getElementsByTagName('specData').item(0).getChildNodes.item(0).getData);
    
       data_char_cell = {};
       while ~isempty(data_char)
           [data_char_cell{end+1},data_char] = strtok(data_char,';');
       end
       
       num_cols = length(strfind(data_char_cell{1},','))+1;
       data_cell = cell(length(data_char_cell),num_cols);
       for i = 1:length(data_char_cell)
           data_char_row = data_char_cell{i};
           data_cell_row = {};
           while ~isempty(data_char_row)
               [data_cell_row{end+1},data_char_row] = strtok(data_char_row,',');
           end
           data_cell(i,:) = deal(data_cell_row);
       end
       
       dc = cellfun(@str2num,data_cell, 'UniformOutput', false);
       mask = cellfun(@isempty,dc);
       data_cell(~mask) = dc(~mask);
       
       fileList_struct.specData = data_cell;
    
    end
       
    
catch
    errordlg('XML file does not contain all fields - stopping!','Import error')
    return
end

try
    spec_node = fileList.getElementsByTagName('spec').item(0);
    instances = spec_node.getElementsByTagName('instance');
    
    for j = 1:instances.getLength;
        fileList_struct.spec(j).rawFile = char(instances.item(j-1).getElementsByTagName('rawFile').item(0).getChildNodes.item(0).getData);
        fileList_struct.spec(j).ID = char(instances.item(j-1).getElementsByTagName('ID').item(0).getChildNodes.item(0).getData);
        if ~isempty(instances.item(j-1).getElementsByTagName('sampleID').item(0).getChildNodes.item(0))
            fileList_struct.spec(j).sampleID = char(instances.item(j-1).getElementsByTagName('sampleID').item(0).getChildNodes.item(0).getData);
        end
    end
    
catch
    errordlg('XML fileList does not contain information about spectral raw files... quitting','Import error')
    return
end

%add some conversions to help this xml format fit with the expected
%fileList structure variable names from legacy code:

%SumTransients conversions
fileList_struct.avTransDir = fileList_struct.avgdTransDirectory;
fileList_struct.setUID = fileList_struct.setIdentifier;

fileList_struct.rawShort = {fileList_struct.spec.ID};
for i=1:length(fileList_struct.rawShort)
    fileList_struct.rawFull{i} = [fileList_struct.rootDirectory,fileList_struct.rawDirectory, fileList_struct.rawShort{i},'.raw'];
end
fileList_struct.rawContig = 0; % set to 1 if raw files all combine to form a single spectrum (never used)

fileList_struct.datDir = [fileList_struct.rootDirectory, fileList_struct.datDirectory];
for i=1:length(fileList_struct.rawShort)
    fileList_struct.datStem{i} = [fileList_struct.rawShort{i},'_'];
end

for i = 1:length(fileList_struct.rawShort)
    fileList_struct.scans{i} = 'all'; %this is a prime example of something that is being kept despite having apparently no function.
end

return
end



function [iOut] = PeakListCluster_0_3(mzIn, maxSpread)
% mzIn = SORTED (ascending) peak m/z values for N spectra
%           (cell length N, each cell is matrix size 1 x n peaks)
% maxSpread = maximum spread of peaks in cluster (m/z or ppm)

% iOut = cluster index of each peak, same size structure as mzIn
%           NB: peak count can be > 1 for one spectrum!  (Means close peaks
%           in spec.)

% MODIFICATION HISTORY
% 
% 14/Nov/07 (v0.1) Added rethrow error to pdist catch
% Rearranged conversion of pdist's to ppm to reduce memory requirements.
% 15/Nov/07 (v0.1b) Added PARAMS.BLOCKSIZE for better splitting up of
% search space.
% 7/Dec/07  (v0.1c) Added flag to speed processing and conserve memory 
% in cases where peak must be present in all spectra
% 17/Dec/07         Removed flag
% 18/Feb/08 (v0.2)  Fixed bug where block boundaries coinciding with peaks
%                   cause false peak clusters
% 19/Feb/08 (v0.3)  Simplified, takes only m/z values and returns cluster
%                   index of each m/z value (and mean cluster m/z)
% 11/Mar/08         Gives hint that PLS toolbox might cause cluster to fail
% 15/May/08         lin 64, bufl and bufh increased from 2* to 5* to try
%                   remove occurrences of failing
% 15/Jul/08         Added hint to rehash toolboxes and path list if cluster
%                   fails (seems function shadowing sometimes needs update)

%% start
MAXBLOCKSIZE = 2000;  % approx. maximum number of peaks in each block during clustering (reduce to solve MEMORY problems)
                      % a good value is 2000

numSpec = length(mzIn);

% mzOut = [];   % output of cluster mean m/z values

% iOut has same size as mzIn
iOut = {};
for i=1:numSpec
    iOut{i} = []; %zeros(size(mzIn{i}));
end

% search over limited space to conserve memory: limit number of peaks in
% block. Find suitable block end points.
mz = unique([mzIn{:}]);   % all data points
i = 1;    % index to mz
bc = 1;   % boundary count
bound = mz(i);    % first mz boundary value
while i < length(mz)
    i = i+MAXBLOCKSIZE;
    if i > length(mz), i=length(mz); end
    bc = bc+1;    
    bound(bc) = mz(i);
end

fprintf('Splitting search space into blocks with max ~%d peaks\n',MAXBLOCKSIZE);
fprintf('Clustering %d input peak lists\n',numSpec);

cc = 0; % cluster counter
tic;
fprintf('Progress:  0%%');
for bc = 1:length(bound)-1
    % cluster m/z values within current boundary
    MZb = [];
    SMb = [];   % spectral membership of each m/z value
    % add buffer to boundaries (to burn later) in case boundary coincides with cluster
    bufl = 5*maxSpread*bound(bc)/1e6; % need 5*maxSpread
    bufh = 5*maxSpread*bound(bc+1)/1e6;
    % find all the peaks from each spectrum within current boundaries
    for i=1:numSpec
        mz = mzIn{i};
        idx = find(bound(bc)-bufl <= mz & mz <= bound(bc+1)+bufh);
        MZb = [MZb mzIn{i}(idx)];
        SMb = [SMb repmat(i,1,length(idx))];
    end
    % catch cases where linkage won't work
    if length(MZb)<3, continue; end
    % cluster
    try
        Y = pdist(MZb.','cityblock');
    catch
        fprintf('Failed on pdist: is Statistics toolbox installed?\n');
        rethrow(lasterror);
    end
    % convert absolute distances to ppm (of mean of pair) distances
    k = 0;
    for j=1:length(MZb)-1
        idx = k+1:k+length(MZb)-j;
        Y(idx) = Y(idx) ./ mean([repmat(MZb(j),1,length(MZb)-j); MZb(j+1:end)]) * 1e6;
        k = k+length(MZb)-j;
    end
    Z = linkage(Y,'complete');  %furthest distance
    clear Y
    % cluster to give the cluster membership of each value in MZb
    try
        CMb = cluster(Z,'criterion','distance','cutoff',maxSpread);
    catch
        fprintf('''cluster'' function failed: is PLS_Toolbox installed?\n');
        fprintf('\tCheck PLS_Toolbox paths are below stats toolbox.\n');
        fprintf('\tIf you still see this error, try:\n');
        fprintf('\t''rehash pathreset''\n');
        fprintf('\t''rehash toolboxreset''\n');
        rethrow(lasterror);
    end
    clear Z
    % calculate mean m/z for each cluster, assuming clusters labelleled
    % monotonically
    CMZb = zeros(1,max(CMb));
    for i=1:max(CMb)
        CMZb(i) = mean(MZb(CMb==i));
    end
    %if cluster m/z is within the boundary buffer, burn it (it will
    %be included in the next region)
    if (bc > 1) && (bc < (length(bound)-1))
        % middle bound
        cBurn = find(bound(bc) > CMZb | CMZb >= bound(bc+1));
    elseif length(bound) == 2
        % single boundary region - do nothing
        cBurn = [];
    elseif bc == length(bound) - 1
        % final bound
        cBurn = find(bound(bc) > CMZb);
    elseif bc == 1
        % first bound
        cBurn = find(CMZb >= bound(bc+1));
    end
    % burn clusters
%     CMZb(cBurn) = 0;
    for i=1:length(cBurn)
        idx = find(CMb==cBurn(i));
        CMb(idx) = [];
        SMb(idx) = [];
        % re-label highest cluster number with cluster just burnt
        % to avoid missing values in the final output
        cl = unique(CMb); % cluster list
        if max(cl) > cBurn(i)
            % cluster membership list
            CMb(CMb==max(cl)) = cBurn(i);
            % and in case in burn list
            cBurn(cBurn==max(cl)) = cBurn(i);
%             % cluster m/z list
%             CMZb(cBurn(i)) = CMZb(max(cl));
%             CMZb(max(cl)) = 0;
        end
    end
%     CMZb = CMZb(CMZb>0);    % remove burned cluster m/z values
    % update output of cluster centers
%     mzOut = [mzOut CMZb];
    % offset cluster numbers so cluster numbers remain unique
    ccOld = cc;
    cc = cc + max(CMb);
    CMb = CMb + ccOld;
    % store cluster membership to output
    for i=1:numSpec
        iOut{i} = [iOut{i} CMb(SMb==i).'];
    end
    if round(toc) || bc==length(bound)-1
        fprintf('\b\b\b');
        fprintf('%2.0f%%',floor(bc/(length(bound)-1)*100));
        tic;
    end
end
fprintf('\n');

cl = unique([iOut{:}]); % cluster list
nc = length(cl);   %number of clusters

% check cluster numbered monotonically
if max(cl) ~= nc
    error('Clustering failed: non-monotonic cluster numbers detected');
end
% check dimensions match
for i=1:numSpec
    if ~isequal(size(mzIn{i}),size(iOut{i})), error(['Clustering failed i=',num2str(i)]); end
end
% % check number of clusters = length of mzBar
% if nc ~= length(find(mzOut))
%     error('Clustering failes');
% end

return
end


function [MZ, MZbar, Y, Ybar, N] = Clusters2Matrices_0_2(Ci, MZc, Yc, minPeaks, fOptions)
% Processes the output from PeakListCluster: expands the cluster labelling
% into matrices MZ, Y, N.
% Also filters the output according to minPeaks and fOptions.
% Sorts the result by MZbar
% 
% INPUT
% Ci - cell vector (length C) of cluster indices for each spectrum
% MZc - cell vector (length C) of spectrum m/z's
% Yc - cell vector (length C) of spectrum y values
% minPeaks - minimum number of peaks for n-filtering
% fOptions - mFilt 1 = remove CLUSTERS containing overlapping peaks =
%                      peaks from same spectrum occurring multiple times
%                  0 = use right-most peak in cluster)
%            nFilt 1 = filter peaks by presence in minimum number of
%                      spectra
%                  0 = no minimum number of spectra
%            remCl 1 = remove clusters filtered out by mFilt of nFilt
%                  0 = keep clusters but set N to 0 for all peaks caught
%                      by filter)
%             sort 1 = sort output by MZbar
%                  0 = don't
% 
% OUTPUT
% MZ - matrix (size C by P) of m/z values for each of P clusters
% MZbar - mean m/z values (length P) for each cluster
% Y - ", y values
% Ybar - ", mean y values (counting zeros as zeros)
% N - ", binary value indicating peak present (0) or absent (1)
% 
% VERSIONS
% 02/Apr/08 v0.1    Initial version
% 03/Apr/08 v0.2    Filtered clusters removed if remCl = 1, added sort as
%                   option

%% quick check
if (length(Ci) ~= length(MZc)) || (length(MZc) ~= length(Yc))
    error('Input not consistent');
end
numSpectra = length(MZc);

%% create cluster matrices
C = unique([Ci{:}]);    % index of cluster numbers to use
MZ = zeros(numSpectra,length(C));
Y = zeros(numSpectra,length(C));
N = zeros(numSpectra,length(C));
for i=1:numSpectra
    % m/z values
    MZ(i,Ci{i}) = MZc{i};
    % intensity values
    Y(i,Ci{i}) = Yc{i};
    % peak present
    N(i,Ci{i}) = 1;
end
% at this point, any overlapping peaks (multiple peaks present in same
% sample in same cluster) will be included but represented by the RIGHT
% most peak (consequence of assigning same value in MZ to multiple values
% in MZc above).
MZbar = MeanNonZero(MZ,1);
Ybar = mean(Y,1);

%% filter multiple peaks out
% if the same replicate has multiple peaks in a cluster, discard that
% cluster and make user aware until a more suitable way to
% differentiate between overlapping peaks can be found (also, currently,
% overlapping peaks will be poorly quantified)
if fOptions.mFilt
    edges = [C(1)-0.5 C+0.5];
    Cr = [];  % clusters to remove
    for i=1:numSpectra
        % count of each cluster
        cc = hist(Ci{i},edges,2);
        % select clusters to remove
        idx = cc>1;
        Cr = [Cr edges(idx)+0.5];
    end
    % found any?
    if ~isempty(Cr)
        warning('Found %d cluster(s) containing overlapping peaks',length(Cr));
        % yes - remove overlapping clusters?
        if fOptions.remCl
            % yes - remove unwanted clusters from index
            C = setdiff(C,Cr);
        else
            % no - just set cluster N's to zero
            N(:,Cr) = 0;
        end
    end
end

%% filter peaks occurring in minimum number of spectra
% filter on all clusters then combine result with any earlier filtering
if fOptions.nFilt
    % find peaks present in less than minimum number of replicates
    idx = sum(N,1)<minPeaks;
    Cr = find(idx);
    % any?
    if ~isempty(Cr)
        fprintf('Found %d clusters(s) containing insufficient peaks\n',length(Cr));
        % yes - remove cluster?
        if fOptions.remCl
            % yes - remove unwanted clusters from index
            C = setdiff(C,Cr);
        else
            % no - just set cluster N's to zero
            N(:,Cr) = 0;
        end
    end
end

%% apply filters to matrices
if fOptions.remCl
    % remove clusters if required
    MZ = MZ(:,C);
    MZbar = MZbar(C);
    Y = Y(:,C);
    Ybar = Ybar(:,C);
    N = N(:,C);
end

%% sort
if fOptions.sort
    % sort by increasing m/z
    [MZbar, idx] = sort(MZbar,'ascend');
    Y = Y(:,idx);
    Ybar = Ybar(:,idx);
    MZ = MZ(:,idx);
    N = N(:,idx);
end

return
end




function y = MeanNonZero(x,dim)
%MEAN   Average or mean value.
% 
% 23/Oct/07 TGP     Mean of only non-zero elements
%   For vectors, MEAN(X) is the mean value of the elements in X. For
%   matrices, MEAN(X) is a row vector containing the mean value of
%   each column.  For N-D arrays, MEAN(X) is the mean value of the
%   elements along the first non-singleton dimension of X.
%
%   MEAN(X,DIM) takes the mean along the dimension DIM of X. 
%
%   Example: If X = [0 1 2
%                    3 4 5]
%
%   then mean(X,1) is [1.5 2.5 3.5] and mean(X,2) is [1
%                                                     4]
%
%   Class support for input X:
%      float: double, single
%
%   See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.
%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.17.4.3 $  $Date: 2005/05/31 16:30:46 $

if nargin==1, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end

  nzCnt = nzCount(x,dim);
  nzCnt(nzCnt==0) = eps;    %avoid div by 0 errors
  y = sum(x)./nzCnt;
else
  nzCnt = nzCount(x,dim);
  nzCnt(nzCnt==0) = eps;
  y = sum(x,dim)./nzCnt;
end

end


function nzCount = nzCount(x,dim);
% returns count of non-zero elements in x in given dimension

x(x~=0) = 1;
nzCount = sum(x,dim);
end













    
