function ReplicateFilter(fileList_in, html_indir, rf_ppm, rf_minpeaks, rf_forcepeak, outfile_peaks, outdir_peaks, outfile_rf)
% 
%   Inputs:
%       1. fileList_in	: the path to a .xml file in FileListManager format
%	2. html_indir	: a full path to a directory containing files output by Stitch. File names should match those in html_infile. 	
%       3. rf_ppm	: the ppm spread for alignment, default 1.5
%       4. rf_minpeaks	: the min number of reps a peak must be in,
%                           default 1                        
%       5. rf_forcepeak	: whether to force the algorithm to find the
%                           nearest peak if it is missing some reps,
%                           default 0 (1 to force)
%	6. outfile_peaks: the full path to an html file to be used for storing output peak list filenames. 
%	7. outdir_peaks	: the full path to a directory in which to store the peak list files (chosen by Galaxy to be linked to outfile_peaks)
%	8. outfile_rf	: full path to a .txt file for storing more output information.  	  
%           
%   Outputs: 
%       creates various files for use in the Galaxy pipeline.
%
%   R.L.Davidson 12/02/2013 
%
%   Version 3.0
%


%% PARAMETERS
if isa(rf_ppm, 'char')
    rf_ppm = str2num(rf_ppm);
end

if isa(rf_minpeaks, 'char')
    rf_minpeaks = str2num(rf_minpeaks);
end

if isa(rf_forcepeak, 'char')
    rf_forcepeak = str2num(rf_forcepeak);
end

batchList = ImportFileListXML(fileList_in);

SFSDir = [html_indir, filesep]; %using one of Galaxy's chosen 'extra files' folder for the output of the Stitch html file.
RFPLDir = [outdir_peaks, filesep]; %using one of Galaxy's chosen 'extra files' folder to hold the output of Replicate Filter

if ~isdeployed
	%% SORT PATH (MOVE STATS TOOLBOX ABOVE PLSTOOLBOX IF NECESSARY)

	original_path = path; %save original path

	% Attempt a basic cluster, using stats toolbox. 
	try
		Y = pdist(rand(1,100)','cityblock');
		Z = linkage(Y,'complete');
		test = cluster(Z,'criterion','distance','cutoff',rf_ppm);

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
		            test = cluster(Z,'criterion','distance','cutoff',rf_ppm);
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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%MAIN LOOP START%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileOutName_html = {}; %this is to capture the output filenames that will then be placed in an HTML file for Galaxy.

gi=0;   %count of number of final (filtered) spectra being produced
NUMREPS = batchList.nReplicates;
for fi=1:NUMREPS:length(batchList.Samples)
    
    %% %%%%%%%%%%HANDLE FILES
    gi = gi+1;
    clear fileList;
    for i=1:NUMREPS
        fileList.name{i} = batchList.Samples(1,fi+i-1).ID;
    end
    numSpectra = NUMREPS;
    
    %make output directories
    %dFold = batchList.filtPkDir;
    %try mkdir(dFold);
    %catch
    %    error(['Cannot make output directory ',dFold]);
    %end
    
    % output files
    groupName = '';
    for i=fi:(fi+NUMREPS-1)
        groupName = [groupName,'_',batchList.Samples(1,i).ID];
    end
    
    
    fileOutName = [RFPLDir,'RFPL_',groupName,'.txt'];
    fileOutName_html{end+1} = ['RFPL_',groupName,'.txt'];
    %paramsFilename = outfile_params;
    rfOutputName = outfile_rf;
    
    try mkdir(RFPLDir);
    catch
        error(['Cannot make processed transient directory ',RFPLDir]);
    end
    
    
    % display parameters
    fprintf('Maximum peak spread = %dppm\n',rf_ppm);
    fprintf('Number of replicates = %d\n',NUMREPS);
    fprintf('Minimum required spectra with matching peaks = %d out of %d\n\n',rf_minpeaks,length(fileList.name));

    %% %%%%%%%%%% READ IN DATA
    peaks = [];
    for i=1:numSpectra
        fName = [SFSDir,'SFS_',fileList.name{i},'.mat'];
        fprintf('Loading file %s...',fName);
        load(fName);
        % only want peak information
        switch SFSParams.Instrument
            case{'ltqft'}
                if SFSParams.C, error('C parameter: unexpected'); end
            case{'qexactive','orbitrap'}
            case{'solarix'}
        end
        peaks(i).mz = f2mz(SFS.peaksf,SFSParams,SFSParams.Instrument);
        peaks(i).y = SFS.peaksd;
        peaks(i).res = SFS.peaksRes;
        peaks(i).snr = SFS.pSNR;
        if isfield(SFS,'peaksFlag')
            peaks(i).flag = SFS.peaksFlag;
        else
            fprintf('\n');
            warning('peaksFlag field not present: setting to 1');
            peaks(i).flag = true(size(peaks(i).y));
        end
        fprintf('done\n');
    end
    clear SFS SFSParams
    
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
    [Ci] = PeakListCluster_0_3(mz, rf_ppm);
    fOptions.mFilt = true;
    fOptions.nFilt = true;
    fOptions.remCl = true;
    fOptions.sort = true;
    [Mclust, MclustBar, Yclust, YclustBar, Nclust] = Clusters2Matrices_0_2(Ci, mz, d, rf_minpeaks, fOptions);
    fprintf('Found %d peaks\n',length(MclustBar));
    
    %% %%%%%%%%%% FILL IN ZEROS IN INTENSITY MATRIX
    % ...from nearest non-clustered peak in each replicate
    if rf_forcepeak
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
                n = length(find(peaks(i).snr(idxr) > 5 & abs(peaks(i).mz(idxr)-mzc(idxc))./peaks(i).mz(idxr)*1e6 > rf_ppm));
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
filtStr = [num2str(rf_minpeaks),'o',num2str(NUMREPS)];
n = length(MclustBar);
eval(['rfOut.n',filtStr,' = ',num2str(n),';']);
TIC = sum(YclustBar);
eval(['rfOut.TIC',filtStr,' = ',num2str(TIC),';']);
if rf_minpeaks > 2
    MCV = std(Yclust,1,1)./YclustBar;
    MCV = median(MCV);
    eval(['rfOut.MCV',filtStr,' = ',num2str(MCV),';']);
end
% all (eg 3 out of 3)
filtIdx = find(sum(Nclust,1)==NUMREPS);
filtStr = [num2str(NUMREPS),'o',num2str(NUMREPS)];
n = length(filtIdx);
eval(['rfOut.n',filtStr,' = ',num2str(n),';']);
TIC = sum(YclustBar(filtIdx));
eval(['rfOut.TIC',filtStr,' = ',num2str(TIC),';']);
if NUMREPS > 2
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
% fig = figure;
% ax = axes;
% hist(ax,dmzppm(dmzppm<50),1000);
% title('Check for high # of very close peaks (false peak discovery) - inc. peak spread to correct');
% ylabel('count');
% xlabel('proximity to closest peak (ppm m/z)');
% %print(fig, outfile_pdf);
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
filtStr1 = [num2str(rf_minpeaks),'o',num2str(NUMREPS)];
filtStr2 = [num2str(NUMREPS),'o',num2str(NUMREPS)];
if isempty(dir(rfOutputName))
    % if file doesn't already exist, create and write headers
    fid = fopen(rfOutputName,'w');
    fprintf(fid,['SAMPLE\tNUM PKS ',filtStr1,'\tTIC ',filtStr1]);
    if rf_minpeaks > 2
        fprintf(fid,['\tMED CV ',filtStr1,' (%%)']);
    end
    fprintf(fid,['\tNUM PKS ',filtStr2,'\tTIC ',filtStr2]);
    if NUMREPS > 2
        fprintf(fid,['\tMED CV ',filtStr2,' (%%)']);
    end
    fprintf(fid,'\r\n');
else
    % if file exists, append
    fid = fopen(rfOutputName,'a');
    fprintf(fid,['SAMPLE\tNUM PKS ',filtStr1,'\tTIC ',filtStr1]);
    if rf_minpeaks > 2
        fprintf(fid,['\tMED CV ',filtStr1,' (%%)']);
    end
    fprintf(fid,['\tNUM PKS ',filtStr2,'\tTIC ',filtStr2]);
    if NUMREPS > 2
        fprintf(fid,['\tMED CV ',filtStr2,' (%%)']);
    end
    fprintf(fid,'\r\n');
end

fprintf(fid,'%s\t%d\t%.3e',...
    groupName,...
    eval(['rfOut.n',filtStr1]),...
    eval(['rfOut.TIC',filtStr1]));
if rf_minpeaks > 2
    fprintf(fid,'\t%.2f',eval(['100*rfOut.MCV',filtStr1]));
end
fprintf(fid,'\t%d\t%.3e',...
    eval(['rfOut.n',filtStr2]),...
    eval(['rfOut.TIC',filtStr2]));
if NUMREPS > 2
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
fprintf(fid,'<title>Replicate Filter Data - Output</title>');
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
fprintf(fid,'');
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
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

fprintf('Finished.\n');

