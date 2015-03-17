function ReplicateFilter_BB(fileList_in, html_infile, html_indir,  rf_ppm, rf_minpeaks, rf_forcepeak, outfile_peaks, outdir_peaks, outfile_rf)
% 
%   Inputs:
%       1. fileList_in	: the path to a .xml file in FileListManager format
%	2. html_infile	: an html file with only <A HREF=""></A> tags, the anchor reference and text content providing each filename as output by Stitch.
%	3. html_indir	: a full path to a directory containing files output by Stitch. File names should match those in html_infile. 	
%       3. rf_ppm	: the ppm spread for alignment, default 1.5
%       4. rf_minpeaks	: the min number of reps a peak must be in,
%                           default 1                        
%       5. rf_forcepeak	: whether to force the algorithm to find the
%                           nearest peak if it is missing some reps,
%                           default 0 (1 to force)
%	6. outfile_peaks: the full path to an html file to be used for storing output peak list filenames. 
%	7. outdir_peaks	: the full path to a directory in which to store the peak list files (chosen by Galaxy to be linked to outfile_peaks)
%	8. outfile_params: full path to a .txt file for storing parameter values (less important when using Galaxy)
%	9. outfile_rf	: full path to a .txt file for storing more output information.  	  
%           
%   Outputs: 
%       creates various files for use in the Galaxy pipeline.
%
%   R.L.Davidson 12/02/2013 
%
%   Version 3.0
%



%% PARAMETERS
    

PARAMS.MAXSPREAD_PPM = rf_ppm;
PARAMS.MINPEAKS = rf_minpeaks;
PARAMS.FORCE_ALLREP_MEAN_INT = rf_forcepeak;


batchList = import_filelist_xml(fileList_in);
batchList.stitchDir = [html_indir, filesep]; %using one of Galaxy's chosen 'extra files' folder for the output of the Stitch html file.
batchList.filtPkDir = [outdir_peaks, filesep]; %using one of Galaxy's chosen 'extra files' folder to hold the output of Replicate Filter

PARAMS.NUMREPS = batchList.numReps;


%% SORT PATH (MOVE STATS TOOLBOX ABOVE PLSTOOLBOX IF NECESSARY)

original_path = path; %save original path

% Attempt a basic cluster, using stats toolbox. 
try
    Y = pdist(rand(1,100)','cityblock');
    Z = linkage(Y,'complete');
    test = cluster(Z,'criterion','distance','cutoff',PARAMS.MAXSPREAD_PPM);

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
                test = cluster(Z,'criterion','distance','cutoff',PARAMS.MAXSPREAD_PPM);
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




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%MAIN LOOP START%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileOutName_html = {}; %this is to capture the output filenames that will then be placed in an HTML file for Galaxy.

gi=0;   %count of number of final (filtered) spectra being produced
for fi=1:PARAMS.NUMREPS:length(batchList.rawShort)

%% %%%%%%%%%%HANDLE FILES
gi = gi+1;
clear fileList;
for i=1:PARAMS.NUMREPS
    fileList.name{i} = batchList.rawShort{fi+i-1};
end
numSpectra = PARAMS.NUMREPS;

%make output directories
%dFold = batchList.filtPkDir;
%try mkdir(dFold);
%catch
%    error(['Cannot make output directory ',dFold]); 
%end

% output files
groupName = '';
for i=fi:(fi+PARAMS.NUMREPS-1)
    groupName = [groupName,'_',batchList.rawShort{i}];
end


fileOutName = [batchList.filtPkDir,'filteredPeaks_',batchList.setUID,groupName,'.txt'];
fileOutName_html{end+1} = ['filteredPeaks_',batchList.setUID,groupName,'.txt'];
%paramsFilename = outfile_params;
rfOutputName = outfile_rf;

try mkdir(batchList.filtPkDir);
catch 
    error(['Cannot make processed transient directory ',batchList.filtPkDir]); 
end


% display parameters
fprintf('Maximum peak spread = %dppm\n',PARAMS.MAXSPREAD_PPM);
fprintf('Number of replicates = %d\n',PARAMS.NUMREPS);
fprintf('Minimum required spectra with matching peaks = %d out of %d\n\n',PARAMS.MINPEAKS,length(fileList.name));

%% %%%%%%%%%% READ IN DATA
peaks = [];
for i=1:numSpectra
    fName = [batchList.stitchDir,'specOut_',batchList.setUID,'_',fileList.name{i},'.mat'];
    fprintf('Loading file %s...',fName);
    load(fName);
    % only want peak information
    if specOutParams.C, error('C parameter: unexpected'); end
    peaks(i).mz = f2mz(specOut.peaksf,specOutParams);
    peaks(i).y = specOut.peaksd;
    peaks(i).res = specOut.peaksRes;
    peaks(i).snr = specOut.pSNR;
    if isfield(specOut,'peaksFlag')
        peaks(i).flag = specOut.peaksFlag;
    else
        fprintf('\n');
        warning('peaksFlag field not present: setting to 1');
        peaks(i).flag = true(size(peaks(i).y));
    end
    fprintf('done\n');
end
clear specOut specOutParams

%% %%%%%%%%%% CLUSTER
fprintf('Clustering peaks...');
d = {};
mz = {};
for i=1:numSpectra
    [mz{i},idx] = sort(peaks(i).mz(peaks(i).flag),'ascend'); % concat'ed m/z values
    % only consider peaks with peaksFlag set when searching for clusters
    d{i} = peaks(i).y(peaks(i).flag);  % concat'ed data values
    d{i} = d{i}(idx);
end
[Ci] = PeakListCluster_0_3(mz, PARAMS.MAXSPREAD_PPM);
fOptions.mFilt = true;
fOptions.nFilt = true;
fOptions.remCl = true;
fOptions.sort = true;
[Mclust, MclustBar, Yclust, YclustBar, Nclust] = Clusters2Matrices_0_2(Ci, mz, d, PARAMS.MINPEAKS, fOptions);
fprintf('Found %d peaks\n',length(MclustBar));

%% %%%%%%%%%% FILL IN ZEROS IN INTENSITY MATRIX
% ...from nearest non-clustered peak in each replicate
if PARAMS.FORCE_ALLREP_MEAN_INT
    for i=1:numSpectra
        % cluster peaks not present in this replicate
        idx = find(Nclust(i,:)==0);
        % found any?
        if ~isempty(idx)
            % what m/z are they in the cluster
            mzc = MclustBar(idx);
            % remove replicate peaks already used in a cluster
            mzr = peaks(i).mz;
            [mzr, idx2] = setdiff(mzr,Mclust(i,:));
            % find closest in the replicate (including noise peaks, regardless of
            % distance to cluster centre)
            [idxr,idxc] = FindClosest_0_3(mzr,mzc,[],'off');
            % set idxr relative to peaks(i)
            idxr = idx2(idxr);
            % warn if high SNR and large distance
            n = length(find(peaks(i).snr(idxr) > 5 & abs(peaks(i).mz(idxr)-mzc(idxc))./peaks(i).mz(idxr)*1e6 > PARAMS.MAXSPREAD_PPM));
            if n, warning('%d non-clustered peaks found with SNR>5 and dmz>MAXPSREAD_PPM',n); end
            % enter into matrix
            Yclust(i,idx(idxc)) = peaks(i).y(idxr);
        end
    end
    % check no more zeros left
    if length(find(Yclust==0)), error('%d zeros still in intensity matrix',length(find(Yclust==0))); end
    % mean intensity for each cluster
end
YclustBar = mean(Yclust,1);

%% %%%%%%%%%% STATISTICS
% output information for user
% filtered (eg 2 out of 3)
% eg: rfOut.n2o3 = length(MclustBar);    % number peaks in at least 2 of 3 spectra
filtStr = [num2str(PARAMS.MINPEAKS),'o',num2str(PARAMS.NUMREPS)];
n = length(MclustBar);
eval(['rfOut.n',filtStr,' = ',num2str(n),';']);
TIC = sum(YclustBar);
eval(['rfOut.TIC',filtStr,' = ',num2str(TIC),';']);
if PARAMS.MINPEAKS > 2
    MCV = std(Yclust,1,1)./YclustBar;
    MCV = median(MCV);
    eval(['rfOut.MCV',filtStr,' = ',num2str(MCV),';']);
end
% all (eg 3 out of 3)
filtIdx = find(sum(Nclust,1)==PARAMS.NUMREPS);
filtStr = [num2str(PARAMS.NUMREPS),'o',num2str(PARAMS.NUMREPS)];
n = length(filtIdx);
eval(['rfOut.n',filtStr,' = ',num2str(n),';']);
TIC = sum(YclustBar(filtIdx));
eval(['rfOut.TIC',filtStr,' = ',num2str(TIC),';']);
if PARAMS.NUMREPS > 2
    MCV = std(Yclust(:,filtIdx),1,1)./YclustBar(filtIdx);
    MCV = median(MCV);
    eval(['rfOut.MCV',filtStr,' = ',num2str(MCV),';']);
end

%% histogram of m/z (ppm) spacings between adjacent peaks
dmz = MclustBar(2:end) - MclustBar(1:end-1);
dmzppm = dmz ./ MclustBar(1:end-1) * 1e6;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This needs to be added - either by exporting to a csv or by allowing 3D acceleration for the VM!
%fig = figure;
%ax = axes;
%hist(ax,dmzppm(dmzppm<50),1000);
%title('Check for high # of very close peaks (false peak discovery) - inc. peak spread to correct');
%ylabel('count');
%xlabel('proximity to closest peak (ppm m/z)');
%print(fig, outfile_pdf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%% OUTPUT
% peak list
fprintf('Saving peaks to %s...',fileOutName);
fid = fopen(fileOutName,'w');
fprintf(fid,'M/Z\tINTENSITY\tNUM SPECTRA (PEAK FLAGGED)\r\n');
for i=1:length(MclustBar)
    fprintf(fid,'%.7f\t%.6e\t%d\r\n',MclustBar(i),YclustBar(i),sum(Nclust(:,i)));
end
fclose(fid);
fprintf('done\n');

% replicate filter output
fprintf('Saving statistics to %s...',rfOutputName);
% filter names
filtStr1 = [num2str(PARAMS.MINPEAKS),'o',num2str(PARAMS.NUMREPS)];
filtStr2 = [num2str(PARAMS.NUMREPS),'o',num2str(PARAMS.NUMREPS)];
if isempty(dir(rfOutputName))
    % if file doesn't already exist, create and write headers
    fid = fopen(rfOutputName,'w');
    fprintf(fid,['SAMPLE\tNUM PKS ',filtStr1,'\tTIC ',filtStr1]);
    if PARAMS.MINPEAKS > 2
        fprintf(fid,['\tMED CV ',filtStr1,' (%%)']);
    end
    fprintf(fid,['\tNUM PKS ',filtStr2,'\tTIC ',filtStr2]);
    if PARAMS.NUMREPS > 2
        fprintf(fid,['\tMED CV ',filtStr2,' (%%)']);
    end
    fprintf(fid,'\r\n');
else
    % if file exists, append
    fid = fopen(rfOutputName,'a');
    fprintf(fid,['SAMPLE\tNUM PKS ',filtStr1,'\tTIC ',filtStr1]);
    if PARAMS.MINPEAKS > 2
        fprintf(fid,['\tMED CV ',filtStr1,' (%%)']);
    end
    fprintf(fid,['\tNUM PKS ',filtStr2,'\tTIC ',filtStr2]);
    if PARAMS.NUMREPS > 2
        fprintf(fid,['\tMED CV ',filtStr2,' (%%)']);
    end
    fprintf(fid,'\r\n');
end

fprintf(fid,'%s\t%d\t%.3e',...
    groupName,...
    eval(['rfOut.n',filtStr1]),...
    eval(['rfOut.TIC',filtStr1]));
if PARAMS.MINPEAKS > 2
    fprintf(fid,'\t%.2f',eval(['100*rfOut.MCV',filtStr1]));
end
fprintf(fid,'\t%d\t%.3e',...
    eval(['rfOut.n',filtStr2]),...
    eval(['rfOut.TIC',filtStr2]));
if PARAMS.NUMREPS > 2
    fprintf(fid,'\t%.2f',eval(['100*rfOut.MCV',filtStr2]));
end
fprintf(fid,'\r\n');
fclose(fid);
fprintf('done\n');

end % batch list loop

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MAIN LOOP END %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%save HTML file for Galaxy

%create html output for Galaxy
fid = fopen(outfile_peaks, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Sum Transient Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>Sum Transient Data - Output</b></h2></div>');
fprintf(fid,'<hr></hr>')
for i=1:length(fileOutName_html)
  	html_code = ['<div style="padding: 3px"><b><a href="',fileOutName_html{i},'">',fileOutName_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'<p>');
fprintf(fid,'Note: ---');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);

%% %%%%%%%%%%%save parameters
%fprintf('Saving parameters to file: %s...\n',paramsFilename);
%fid = fopen(paramsFilename,'w');
%if ~fid, error('Cannot create message file'); end
%fprintf(fid,'FILE VERSION:\t%s',mfilename);
%fprintf(fid,'\r\n\r\nCLUSTER PARAMETERS:\r\n');
%fn = fieldnames(PARAMS);
%for i=1:length(fn)
%    val = PARAMS.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else
%        fc=fprintf(fid,'\r\n\t%s:',fn{i});
%        fc=fc & fprintf(fid,'\t%.10g',val);
%    end
%    if ~fc, error('Cannot write to message file'); end
%end
%fclose(fid);

%% RETURN PATH TO ORIGINAL STATE (In case of rejigging due to cluster function issues)

path(original_path);
rehash pathreset;
rehash toolbox;



fprintf('Finished.\n');


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


function mz = f2mz(f, P, y)
%converts from freq, returns m/z
%P.A, P.B and P.C contain calibration parameters
%if C is supplied, it is expected to apply to the form of calibration equation derived by Masselon et al. 2002, ref 63:
% m/z = A./f + B./f.^2 + C*y./f.^2, where y is the intensity of each peak

numTerms = 2;   % number of parameters in calibration equation

if nargin == 1
    % no P parameter supplied
    disp('WARNING: Insufficient parameters specified, using nominal values');
    numTerms = 2;
    A = 107368.5e3;
    B = -750.461e6;
elseif isempty(P)
    % P is empty: use defaults without warning
    numTerms = 2;
    A = 107368.5e3;
    B = -750.461e6;
elseif ~isfield(P,'C') || P.C==0
    % no or zero C (therefore y irrelevant) parameter: use two-term calibration
    numTerms = 2;
    A = P.A;
    B = P.B;
    C = 0;
else
    % 3-term calibration
    numTerms = 3;
    A = P.A;
    B = P.B;
    C = P.C;
end

switch numTerms
    case 2
        mz = A./f + B./f.^2;
    case 3
        mz = A./f + B./f.^2 + C*y./f.^2;
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

function nzCount = nzCount(x,dim);
% returns count of non-zero elements in x in given dimension

x(x~=0) = 1;
nzCount = sum(x,dim);



function [idxTarget,idxIn] = FindClosest(target,inputList,maxErr,errType)
%find THE (if present) closest target value to EACH inputList value
%ie. length(idxTarget) = length(idxIn) = length(inputList)
%returns the indices (if found) of the points in the inputList and the
%target list
%maxErr is maximum distance allowed from points
%errType = 'abs' or 'ppm' or 'off' for maxErr given as absolute or ppm or
%no maxErr.

%Updates:
%09/05/07, version 0.1: lists can be unsorted
%26/May/07, v0.1:   update for length(target)=1 and length(inputList)=1
%07/Dec/07, v0.2:   now returns only THE closest value, provided it is within
%                   the maxErr distance. Previously returned ALL closest values!
%18/Feb/08  v0.3:   minor bug correction in event only single target value
%                   (NB v0.2 seems to not have been renamed)

%check for empty input
if (isempty(inputList) | isempty(target)); warning('Empty input'); idxTarget=[]; idxIn=[]; return; end;

%sort inputList
[inputList,idxSortInput] = sort(inputList.');
inputList = inputList.';

%sort target
[target,idxSortTarget] = sort(target.');
target = target.';

%find index of closest in target list to each input
if length(target)==1
    ti = ones(size(inputList)); % each target is closest to input
else
    % interpolate inputList onto target
    ti = interp1(target, 1:length(target), inputList, 'nearest', 'extrap');
end

%find those with distance error <= maxErr
switch errType
    case 'off'
        ii = 1:length(inputList);
    case 'abs'
        e = abs(inputList - target(ti));
        ii = find(e <= maxErr);
        ti = ti(ii);
    case 'ppm'
        e = abs(inputList - target(ti)) ./ target(ti);
        ii = find(e <= maxErr/1e6);
        ti = ti(ii);
    otherwise
        error('Format choise not recognised');
end

%return indices wrt original (unsorted) input
idxTarget = idxSortTarget(ti).';
idxIn = idxSortInput(ii).';
