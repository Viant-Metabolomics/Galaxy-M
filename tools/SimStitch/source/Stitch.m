function Stitch(fileList, html_indir, ListofCalibrants, csv_str, html_outfile_stitch, html_outdir_stitch, html_outfile_peaks, html_outdir_peaks, noise_file, align_file)
%
%   Stitch(config_file)
%
%       inputs (the config file should be a .csv file, with each variable on a separate row and an identifying number in the first column:
%           1,     SET = set number from FileListManager
%           
%           2,     sOptions.files.userawonly_on
%           3,     sOptions.noise.noisefilt_on
%           4,     sOptions.noise.includenoisepeaks_on
%           5,     sOptions.noise.min_snr
%           6,     sOptions.noise.noiseremoveallknown_on
%           7,     sOptions.noise.highnoiseknownregions_mz
%           8,     sOptions.cal.on
%           9,     sOptions.cal.blank_on
%           10,    sOptions.cal.mode
%           11,    sOptions.cal.weighted_on
%           12,    sOptions.cal.min_range
%           13,    sOptions.cal.max_pk_d
%           14,    sOptions.cal.min_snr
%           15,    sOptions.cal.ac_m
%           16,    sOptions.cal.ac_c
%           17,    sOptions.cal.ac_p
%           18,    sOptions.align.ms_align_on
%           19,    sOptions.align.intensity_correct_on
%           20,    sOptions.align.m_min_range
%           21,    sOptions.align.min_snr
%           22,    sOptions.align.max_pk_d
%           23,    sOptions.segrange

%
%       outputs:
%           outputs various files to stitchDir specified in FileListManager
%
%
%   R.L. Davidson 14/02/2013 updated from original T.S.Payne
%
%   Version 3.0


% ########################################################################
%%  PARAMETERS
% ########################################################################

%% 1. Read options from csv file csv_str

tmpName = tempname;
fid = fopen(tmpName,'w');
fprintf(fid, csv_str);
fclose(fid);

try 
    configs = csvread(tmpName);
catch err
    if strcmp(err.identifier,'MATLAB:textscan:BadFormatString')
        sprintf('Stitch config file may be empty!')
    elseif strcmp(err.identifier,'MATLAB:textscan:handleErrorAndShowInfo')
        sprintf('Stitch config file may contain text. Please use only numbers separated by commas.')
    else
        sprintf('Unidentified problem reading config file. Err code: %s',err.identifier)
    end


	%make an output for Galaxy
	fid = fopen(html_outfile_peaks,'w');
	fprintf(fid,'Stitch.m failed at reading config.csv file!');
	fclose(fid);

	fid = fopen(html_outfile_stitch,'w');
	fprintf(fid,'Stitch.m failed at reading config.csv file!');
	fclose(fid);

    return
end

config_ind = configs(:,1); %the first column in configs.csv should be numbers that identify the variable in subsequent columns.

%% 2. %% Validating and assigning the full list of possible arguments with
% reference to config_file -> configs


temp = find(config_ind==2); %USERAWONLY_ON ######################################   2
if ~isempty(temp)
    FILES.USERAWONLY_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    FILES.USERAWONLY_ON = 0; %DEFAULT
    % use only the raw files (ie not processed transient data) in stitching process?
end
if length(FILES.USERAWONLY_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 2.Too many values.')
    return
elseif isempty(FILES.USERAWONLY_ON)
    FILES.USERAWONLY_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end

temp = find(config_ind==3); %NOISE.NOISEFILT_ON #################################   3
if ~isempty(temp)
    NOISE.NOISEFILT_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    NOISE.NOISEFILT_ON = 1; %DEFAULT
    %measure and filter noise?
end
if length(NOISE.NOISEFILT_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 3.Too many values.')
    return
elseif isempty(NOISE.NOISEFILT_ON)
    NOISE.NOISEFILT_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end



temp = find(config_ind==4); %NOISE.INCLUDENOISEPEAKS_ON   ######################    4
if ~isempty(temp)
    NOISE.INCLUDENOISEPEAKS_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    NOISE.INCLUDENOISEPEAKS_ON = 1; %DEFAULT
    %include noise peaks (peaks with data points below data threshold) in output and flag them?
end
if length(NOISE.INCLUDENOISEPEAKS_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 4.Too many values.')
    return
elseif isempty(NOISE.INCLUDENOISEPEAKS_ON)
    NOISE.INCLUDENOISEPEAKS_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end



temp = find(config_ind==5); %NOISE.MIN_SNR  ####################################    5    
if ~isempty(temp)
    NOISE.MIN_SNR = configs(temp,2:max(find(configs(temp,1:end))));
else    
    NOISE.MIN_SNR = 10; %DEFAULT
    %peak must have a height > NOISE.MIN_SNR*(standard deviation of complex noise data points).
    %setting MIN_SNR to 3.5 roughly equates to previous NOISE.NFMULT of 2.5
end
if length(NOISE.MIN_SNR)>1   %VALIDATE 1
    sprintf('Error: config variable number 5.Too many values.')
    return
elseif isempty(NOISE.MIN_SNR)
    NOISE.MIN_SNR=0; %The line above ignores zeros, but was useful for checking extra values. 
end



temp = find(config_ind==6); %NOISE.NOISEREMOVEALLKNOWN_ON   ####################    6
if ~isempty(temp)
    NOISE.NOISEREMOVEALLKNOWN_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    NOISE.NOISEREMOVEALLKNOWN_ON = 1; %DEFAULT
    %set to 1 to remove all known highly noisy regions
end
if length(NOISE.NOISEREMOVEALLKNOWN_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 6.Too many values.')
    return
elseif isempty(NOISE.NOISEREMOVEALLKNOWN_ON)
    NOISE.NOISEREMOVEALLKNOWN_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end



temp = find(config_ind==7); %HIGHNOISEKNOWNREGIONS_MZ   ########################    7
if ~isempty(temp)
    HIGHNOISEKNOWNREGIONS_MZ = configs(temp,2:max(find(configs(temp,1:end))));
else    
    HIGHNOISEKNOWNREGIONS_MZ = []; %DEFAULT
    %known regions of high noise (check spectra)
end
if mod(length(HIGHNOISEKNOWNREGIONS_MZ),2)   %VALIDATE 1
    sprintf('Error: config variable number 7. Values must be paired.')
    return
else
    HIGHNOISEKNOWNREGIONS_MZ = reshape(HIGHNOISEKNOWNREGIONS_MZ,2,length(HIGHNOISEKNOWNREGIONS_MZ)/2)'; % pair up high noise regions
end  



temp = find(config_ind==8); %CAL.ON     ########################################    8
if ~isempty(temp)
    CAL.ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.ON = 1; %DEFAULT
    %internally recalibrate (1) or not (0) where possible
end
if length(CAL.ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 8. Too many values.')
    return
elseif isempty(CAL.ON)
    CAL.ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end 



temp = find(config_ind==9); %CAL.BLANK_ON   ####################################    9
if ~isempty(temp)
    CAL.BLANK_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.BLANK_ON = 1; %DEFAULT
    %internally recalibrate blank spectra (1) or use external calibration only (0): has no effect if CAL.ON = 0
end
if length(CAL.BLANK_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 9. Too many values.')
    return
elseif isempty(CAL.BLANK_ON)
    CAL.BLANK_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end 


temp = find(config_ind==10); %CAL.MODE      #####################################   10
if ~isempty(temp)
    CAL.MODE = configs(temp,2:max(find(configs(temp,1:end))));
    if length(CAL.MODE)>1   %VALIDATE 1
        sprintf('Error: config variable number 10. Too many values.')
        return
    elseif isempty(CAL.MODE)
        CAL.MODE=0; %The line above ignores zeros, but was useful for checking extra values. 
    end 
    
    switch CAL.MODE         %VALIDATE 2
        case 1
            CAL.MODE = 'ab';
        case 2
            CAL.MODE = 'df';
        case 3
            CAL.MODE = 'abc';
        case 4
            CAL.MODE = 'linear';
        otherwise
            sprintf('Error: config variable number 10 should be in range 1:4. 1=ab,2=df,3=abc,4=linear.')
            return
    end
else
    CAL.MODE = 'ab'; %DEFAULT
    %method of calibration ('df','ab','abc','linear') (freq. shifting, ab(c) parameter, linear optimisation)
end  
    
temp = find(config_ind==11); %CAL.WEIGHTED_ON   #################################   11
if ~isempty(temp)
    CAL.WEIGHTED_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.WEIGHTED_ON = 0 %DEFAULT
    %is the calibration best fit weighted to the abundances?
end
if length(CAL.WEIGHTED_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 11. Too many values.')
    return
elseif isempty(CAL.WEIGHTED_ON)
    CAL.WEIGHTED_ON = 0; %The line above ignores zeros, but was useful for checking extra values. 
end   
    

temp = find(config_ind==12); %CAL.MIN_RANGE     ################################    12
if ~isempty(temp)
    CAL.MIN_RANGE = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.MIN_RANGE = 50; %DEFAULT
    %minimum separation (in PERCENT of total range) between any 2 calibration points for A and B parameter calibration, otherwise only B parameter calibrated
end
if length(CAL.MIN_RANGE)>1   %VALIDATE 1
    sprintf('Error: config variable number 12. Too many values.')
    return
elseif isempty(CAL.MIN_RANGE)
    CAL.MIN_RANGE=0; %The line above ignores zeros, but was useful for checking extra values. 
    
elseif CAL.MIN_RANGE>100 || CAL.MIN_RANGE<0 %VALIDATE 2
    sprintf('Error: config variable number 12. Needs to be a percentage value in range 0:100.')
    return
end  


temp = find(config_ind==13); %CAL.MAX_PK_D      ################################    13
if ~isempty(temp)
    CAL.MAX_PK_D = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.MAX_PK_D = 6.5; %DEFAULT
    %maximum distance between measured (externally calibrated) peaks and calibrants
end
if length(CAL.MAX_PK_D)>1   %VALIDATE 1
    sprintf('Error: config variable number 13. Too many values.')
    return
elseif isempty(CAL.MAX_PK_D)
    CAL.MAX_PK_D=0; %The line above ignores zeros, but was useful for checking extra values. 
end      
    

temp = find(config_ind==14); %CAL.MIN_SNR       ################################    14
if ~isempty(temp)
    CAL.MIN_SNR = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.MIN_SNR = 10; %DEFAULT
    %calibration only takes place on peaks with at least this SNR
end
if length(CAL.MIN_SNR)>1   %VALIDATE 1
    sprintf('Error: config variable number 14. Too many values.')
    return
elseif isempty(CAL.MIN_SNR)
    CAL.MIN_SNR=0; %The line above ignores zeros, but was useful for checking extra values. 
end     


temp = find(config_ind==15); %CAL.AC_M      ####################################    15
if ~isempty(temp)
    CAL.AC_M = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.AC_M = 0.05e-3; %DEFAULT
    %window abundance correction 'm' (gradient) parameter
end
if length(CAL.AC_M)>1   %VALIDATE 1
    sprintf('Error: config variable number 15. Too many values.')
    return
elseif isempty(CAL.AC_M)
    CAL.AC_M=0; %The line above ignores zeros, but was useful for checking extra values. 
end 
    
temp = find(config_ind==16); %CAL.AC_C      ####################################    16
if ~isempty(temp)
    CAL.AC_C = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.AC_C = 2e-3;  %DEFAULT
    %window abundance correction 'c' (intercept) parameter
end
if length(CAL.AC_C)>1   %VALIDATE 1
    sprintf('Error: config variable number 16. Too many values.')
    return
elseif isempty(CAL.AC_C)
    CAL.AC_C=0; %The line above ignores zeros, but was useful for checking extra values. 
end 

temp = find(config_ind==17); %CAL.AC_P      ####################################    17
if ~isempty(temp)
    CAL.AC_P = configs(temp,2:max(find(configs(temp,1:end))));
else    
    CAL.AC_P = 0; %DEFAULT
    %window abundance correction pivot point (relative to window start m/z) in m/z
end
if length(CAL.AC_P)>1   %VALIDATE 1
    sprintf('Error: config variable number 17. Too many values.')
    return
elseif isempty(CAL.AC_P)
    CAL.AC_P=0; %The line above ignores zeros, but was useful for checking extra values. 
end
 

temp = find(config_ind==18); %ALIGN.MZ_ALIGN_ON     ############################    18
if ~isempty(temp)
    ALIGN.MZ_ALIGN_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    ALIGN.MZ_ALIGN_ON = 1; %DEFAULT
    %horizontal alignment
end
if length(ALIGN.MZ_ALIGN_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 18. Too many values.')
    return
elseif isempty(ALIGN.MZ_ALIGN_ON)
    ALIGN.MZ_ALIGN_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end


temp = find(config_ind==19); %ALIGN.INTENSITY_CORRECT_ON    ####################    19 
if ~isempty(temp)
    ALIGN.INTENSITY_CORRECT_ON = configs(temp,2:max(find(configs(temp,1:end))));
else    
    ALIGN.INTENSITY_CORRECT_ON = 1; %DEFAULT
    %vertical correction
end
if length(ALIGN.INTENSITY_CORRECT_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 19. Too many values.')
    return
elseif isempty(ALIGN.INTENSITY_CORRECT_ON)
    ALIGN.INTENSITY_CORRECT_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
end

temp = find(config_ind==20); %ALIGN.M_MIN_RANGE     ############################    20
if ~isempty(temp)
    ALIGN.M_MIN_RANGE = configs(temp,2:max(find(configs(temp,1:end))));
else    
    ALIGN.M_MIN_RANGE = 20; %DEFAULT
    %minimum separation (in PERCENT of total range) between any 2 alignment points to use m scaling alignment, otherwise only use c shifting parameter
end
if length(ALIGN.M_MIN_RANGE)>1   %VALIDATE 1
    sprintf('Error: config variable number 20. Too many values.')
    return
elseif isempty(ALIGN.M_MIN_RANGE)
    ALIGN.M_MIN_RANGE=0; %The line above ignores zeros, but was useful for checking extra values. 
    
elseif ALIGN.M_MIN_RANGE>100 || ALIGN.M_MIN_RANGE<0 %VALIDATE 2
    sprintf('Error: config variable number 20. Needs to be a percentage value in range 0:100.')
    return
end


temp = find(config_ind==21); %ALIGN.MIN_SNR     #################################   21
if ~isempty(temp)
    ALIGN.MIN_SNR = configs(temp,2:max(find(configs(temp,1:end))));
else    
    ALIGN.MIN_SNR = 6.5; %DEFAULT
    %alignment only takes place on peaks with at least this SNR
end
if length(ALIGN.MIN_SNR)>1   %VALIDATE 1
    sprintf('Error: config variable number 21. Too many values.')
    return
elseif isempty(ALIGN.MIN_SNR)
    ALIGN.MIN_SNR=0; %The line above ignores zeros, but was useful for checking extra values. 
end

temp = find(config_ind==22); %ALIGN.MAX_PK_D    ################################    22
if ~isempty(temp)
    ALIGN.MAX_PK_D = configs(temp,2:max(find(configs(temp,1:end))));
else    
    ALIGN.MAX_PK_D = 1.5; %DEFAULT
    %maximum span of peaks to use in alignment
end
if length(ALIGN.MAX_PK_D)>1   %VALIDATE 1
    sprintf('Error: config variable number 22. Too many values.')
    return
elseif isempty(ALIGN.MAX_PK_D)
    ALIGN.MAX_PK_D=0; %The line above ignores zeros, but was useful for checking extra values. 
end  

% #####################ADDED FOR MSMS PROCESSING - WINDOW SELECT
temp = find(config_ind==23); %SEG.RANGE    ########################    28
if ~isempty(temp)
    SEG.RANGE = configs(temp,2:max(find(configs(temp,1:end))));
else    
    SEG.RANGE = []; %DEFAULT
    %choose all windows for stitching (normal). 
end
if length(SEG.RANGE)>2   %VALIDATE 1
    sprintf('Error: config variable number 28. Too many values.')
    return
elseif length(SEG.RANGE)==1
    sprintf('Error: config variable number 28. Too few values. Requires variable ID (28) and start and end of range e.g. 28,3,7 or 28,1,6. For full range do not include variable ID 28 or use 28,0,0.') 
    return
elseif isempty(SEG.RANGE)
    SEG.RANGE=[]; %The line above ignores zeros, but was useful for checking extra values. If the user put in 0,0 then this is the value for choosing all windows
end
    
if ~isempty(SEG.RANGE)  %VALIDATE 2  - make sure no one starts the range with 0. e.g 28,0,7 should be 28,1,7
    if SEG.RANGE(1)==0
        SEG.RANGE(1)=1;
    end
    if SEG.RANGE(1)>SEG.RANGE(2) %VALIDATE 3 - make sure that the pair are ordered start and end. i.e. increasing values
        sprintf('Error: config variable number 28. Range needs to have smallest value first. Current range has larger value at start than at end.')
        return
    end
end

% List of calibrants
CAL.FILE = ListofCalibrants;

% ---DEBUG--- %
DEBUG.QUICKPEAKDETECT = 0;    % speed up peak-peaking for debug only (not accurate)
DEBUG.SAVELOADPEAKS = 0;      % save or load intermediate list of detected peaks in _peaks.mat file?

DISPLAY_HDRS = 0;   %display file header info
ZFILLS = 1;        %the number of zero-fills to be applied to transients
X_ZFILLS = 1;      %zero-fills assumed to have been applied by xCalibur for .raw file spectra

% ---FILE--- %
BLANKSTR = 'blank';                 % group name for blank samples

% ---STITCHING--- %
STITCH.ON = 1;

% ---TIC VARS ---%
TICfilter = [];
avTICs= [];
RAWname = {};

%% 3. Setup stitch; get list of files and make folders
%extract the fileList information from the xml file
fileListpath = fileList;
fileList = ImportFileListXML(fileList);

%adjust fileList directories to match Galaxy temporary html file paths
SFSDir = [html_outdir_stitch,filesep]; %this updates the fileList directory instructions for files to be output by Stitch (Galaxy output no 1)
SPLDir = [html_outdir_peaks,filesep]; %this updates the fileList directory instructions for peak files to be output by Stitch (Galaxy output no 2)

%make output directories
dFold = SFSDir; % Stitched mass spectra
[sdir, mess, messid] = mkdir(dFold);
if ~sdir, error(['Cannot make stitched spectra directory ',dFold,': ',messid]); end

dFold = SPLDir; % Stitched peak list
[sdir, mess, messid] = mkdir(dFold);
if ~sdir, error(['Cannot make peak list directory ',dFold,': ',messid]); end

% debug on?
if DEBUG.QUICKPEAKDETECT
    % warn user
    warning('DEBUG.QUICKPEAKDETECT is on: accuracy reduced');
    reply = input('continue - y/n? [n] ','s');
    if isempty(reply) || strcmp(reply,'n')
        error('Stopped by user');
    end
end

%% 4. Main stitching loop over samples

stitchfile_html = {}; %these are place holders to catch the names of output files to be put into ...
peaklist_html = {};  %...html files for Galaxy to use. See collation at end of main loop. 

NOISE.ESTIMATE = 1; % Manually estimate noise from data (default option for data where we can access the transients).

for fi=1:fileList.nDataFiles % loop over each sample
    
    %% a. Setup output + load lable information sample
    % 1. Setup output
    opfname.stitchm = ['SFS_',fileList.Samples(1,fi).ID,'.mat'];
    stitchfile_html{end+1} = opfname.stitchm;
    peaklist_html{end+1} = ['SPL_',fileList.Samples(1,fi).ID,'.txt'];
    %if we want to save the result, check if already stitched    
    fprintf('\nFile: %s\n',fileList.Samples(1,fi).ID); 
    fprintf('Looking for stitched file...');
    fileName = [SFSDir,opfname.stitchm]; %Check stitched file or peak list?
    fileName = dir(fileName);
    if ~isempty(fileName), fprintf('stitched file exists, skipping\n'); continue;
    else fprintf('not found\n'); end
    
    % label whether or not this is a "blank" sample
    if strncmpi(fileList.Samples(fi).sampleID,BLANKSTR,length(BLANKSTR))
        fileList.isBlank(fi) = true;
    else
        fileList.isBlank(fi) = false;
    end
    
    % label whether or not this sample should be internally recalibrated
    % (cal is on and this isn't a blank or it is and we want to
    % calibrate blanks)
    if CAL.ON && (~fileList.isBlank(fi) || (fileList.isBlank(fi) && CAL.BLANK_ON))
        fileList.calOn(fi) = true;
    else
        fileList.calOn(fi) = false;
    end
    
    % 2. Load (uncalibrated) frequency spectra 
    if FILES.USERAWONLY_ON %Thermo QExactive
        % use reduced profile data from .raw files only
        options.getSpec = 1;
        options.sumScans = 1;
        fprintf('Loading reduced profile data from %s.raw\n',fileList.Samples(1,fi).ID);
        % get spectral data and parameters from raw file
        [FS,FSParams,Instrument,NULL_REGION] = GetRawProfileFS_v3(fileList,options,fi,DISPLAY_HDRS,X_ZFILLS);
        NOISE.ESTIMATE = 0; % Use noise values obtained from RAW file
        fileList.Instrument = Instrument;
    else %Thermo FT-Ultra or Bruker Solarix
        fileName = fullfile(html_indir,['FS_',fileList.Samples(1,fi).ID,'.mat']);
        fprintf('Loading pre-processed frequency spectra: %s...',fileName);
        try
            load(fileName);
        catch
            fprintf('not found!\n');
            rethrow(lasterror);
        end
        fprintf('done\n');
        
        % Extract information regarding window overlap from output
        % SumTransients and ProcessTransients
        NULL_REGION = FSParams(1,1).NULL_REGION;
        
        % Determine instrument type
        if strcmp(fileList.Samples(1,1).dataFile(end-2:end),'raw') || strcmp(fileList.Samples(1,1).dataFile(end-2:end),'RAW')
            fileList.Instrument = 'ltqft';
        else
            fileList.Instrument = 'solarix';
        end
    end
        
    % ##########################ADDED FOR MSMS PROCESSING
    if ~isempty(SEG.RANGE) %reduce P and SD to just the chosen range
        if SEG.RANGE(2)>length(P) %check that range is within bounds.
            sprintf('Error: SEG RANGE (2) is greater than total number of windows!')
            return
        else
            FSParams = FSParams(SEG.RANGE(1):SEG.RANGE(2));
            FS = FS(SEG.RANGE(1):SEG.RANGE(2));
        end
    end
    
%     % 18-08-2015 JE: TEMPORARY CODE FOR QE NATURE PROT EXPERIMENTS
%     %delete_windows = 2:length(FS); % Keep first window;
%     delete_windows = 1; % Remove first window
%     FS(delete_windows) = []; FSParams(delete_windows) = []; 
    STITCH.NULL_REGION = NULL_REGION;    
    numWindows = length(FSParams);
    
    %% b. Check data quality
    for si=1:numWindows % Are frequency and intensity vectors of same length?
        if length(FS(si).data)~=length(FS(si).f), error('Spectral intensity and f size mismatch'); end
    end
    
    fprintf('Checking parameters...'); %Is trapping voltage constant
    if strcmp(fileList.Instrument,'ltqft')
        % Trapping voltage always 1 
        if ~FILES.USERAWONLY_ON
            if ~isempty(find([FSParams.VT]~=1, 1))
                error('Unexpected changing or non-unity VT value');
            end
        end
    else
        %ignore VT parameter as it is not produced in the QExactive RAW
        %file or the Bruker .d file for some reason.
        [FSParams.VT] = deal(1); % create VT and set to 1 for each spec.
    end
    
    % Are the calibration parameters OK?
    switch(fileList.Instrument)
        case{'ltqft'}
            if ~isempty(find([FSParams.C]~=0, 1))
                error('Unexpected non-zero C parameter in input files');
            end
            fprintf('done\n');
        case{'qexactive','orbitrap'}
            if ~isempty(find([FSParams.A]~=0, 1))
                error('Unexpected non-zero A parameter in input files');
            end
            fprintf('done\n');
        case{'solarix'}
            %Something for Bruker data
    end
    
    % OLD: sampling frequency correction for ltqft; according to Thermo
    % representatives the peak intensity is not related to the sampling
    % frequency. Most likely this code will be removed at a later stage and
    % be replaced with a normalization step.
    if strcmp(fileList.Instrument,'ltqft')
        if ~FILES.USERAWONLY_ON
            fprintf('Correcting for varying sampling frequency...');
            % find maximum sampling frequency
            fs_max = max([FSParams.BW]) * 2; %BW = bandwidth which is related to the sampling frequency
            
            % correct each SIM window
            for si=1:numWindows
                % gamma is multiplier to achieve maximum fs
                gamma = fs_max / (FSParams(si).BW * 2);
                FS(si).data = FS(si).data * gamma;
            end
            fprintf('done\n');
        end
    end
      
    %% c.  Estimate noise in each window
    if NOISE.ESTIMATE
        %d.1. fit Rayleigh distribution to signal-free region. This region
        %is automatically determined from the data
        fprintf('Measuring noise\n');
        
        for si = 1:numWindows
            % regions of spectrum (in Hz) to NOT use for noise measurement
            exRegion = mz2f(HIGHNOISEKNOWNREGIONS_MZ,FSParams(si),fileList.Instrument); %Update for other machines
            % get data rms noise, sigma of Rayleigh fit, and freq range over which it is measured
            [~,nSigma,nfRange, numPoints] = GetNoiseLevel_0_5(FS(si).data,FS(si).f,exRegion,FSParams(si),fileList.Instrument);
            FSParams(si).dNrange = f2mz(fliplr(nfRange),FSParams(si),fileList.Instrument); %Update for other machines
            % standard deviation of noise (in real, imaginary components of FT,
            % assumed equal and zero-mean Gaussian, estimated from Rayleigh
            % distribution fit of pure noise in signal-free region of
            % amplitude-mode spectrum).
            FSParams(si).dNsd = nSigma;
            % peak height threshold
            FSParams(si).dNthresh = NOISE.MIN_SNR * FSParams(si).dNsd;
            fprintf('dNthresh=%.2f (%.2f x noise sd, range %.3f-%.2f)\n', ...
                FSParams(si).dNthresh,NOISE.MIN_SNR,FSParams(si).dNrange(1),FSParams(si).dNrange(2));
            FSParams(si).dNnumpoints = numPoints; %added by RLD 14/02/2013
        end
    else
        % Use noise estimate from RAW files (ie stitching reduced profile spectra)
        for si = 1:numWindows
            % set irrelevant parameters that are usually output by GetNoiseLevel_0_5
            % to zero
            
            FSParams(si).dNnumpoints = 0;
            FSParams(si).dNrange = [0 0];
        end
        % Filter noise peaks based on SNR values found by getrawprofile
    end
    
    %% d. Threshold data based on estimated noise level; unless noise peaks
    % should be included in output.
    if NOISE.NOISEFILT_ON && ~NOISE.INCLUDENOISEPEAKS_ON
        % apply hard threshold to zero all data below noise; data adjacent
        % to maxima abnove noise threshold are retained because these are
        % used to estimate the peak area
        for si = 1:numWindows
            % Quick estimation of local maxima
            peaksn = maxima(FS(si).data,length(FS(si).data)) + 1; % +1 because output maxima is always 1 value below correct maximum
            % only interested in peaks with maximum data point above the
            % noise threshold
            if NOISE.ESTIMATE
                peaksn = peaksn(FS(si).data(peaksn) > FSParams(si).dNthresh); %dNthresh contains estimated noise threshold 
            else
                peaksn = peaksn(FS(si).SNR(peaksn) > NOISE.MIN_SNR); % Noise filter based on SNR values from getrawfile (for QE)
            end
            fprintf('Window %d, number peaks > thresh = %d\n',si,length(peaksn));
            % and the points either side
            peaksn = unique([peaksn peaksn-1 peaksn+1]);
            % all points above noise threshold
            if NOISE.ESTIMATE
                idx = find((FS(si).data - FSParams(si).dNthresh) >= 0);
            else
                idx = find(FS(si).SNR >= NOISE.MIN_SNR);
            end
            % above noise threshold OR adjacent to peak
            idx = union(idx, peaksn);
            % all other points to zero
            idx = setdiff(1:length(FS(si).data), idx);
            FS(si).data(idx) = 0;
        end
    end
    
    %% e. Detect, measure and flag peaks
        % check for existing _peaks file (debug, saves time to dump peaks)
    fName = [SPLDir,'PL_',fileList.Samples(1,fi).ID,'.mat']; %Path is wrong since peak hasnt been stitched yet.
    files = dir(fName);
    if DEBUG.SAVELOADPEAKS && ~isempty(files)
        fprintf('Loading peak list in file %s\n',fName);
        load(fName);
    else
        if DEBUG.QUICKPEAKDETECT
            peakMethod = 'centroid';
        else
            peakMethod = 'KCEinterp';
        end
        fprintf('Finding peaks:\n');
        for si = 1:numWindows
            fprintf('SIM window %d...',si);
            if NOISE.INCLUDENOISEPEAKS_ON
                if NOISE.ESTIMATE
                    % pick all peaks, including all noise peaks
                    [PL(si).y, PL(si).h, PL(si).f, PL(si).res] = PeakDetectArea_0_3(FS(si).data,FS(si).f,peakMethod,0,FSParams(si),fileList.Instrument,[],[],NOISE); 
                    % flag peaks with height > threshold
                    PL(si).flag = PL(si).h > FSParams(si).dNthresh;
                else
                    % pick all peaks, including all noise peaks
                    [PL(si).y, PL(si).h, PL(si).f, PL(si).res] = PeakDetectArea_0_3(FS(si).data,FS(si).f,peakMethod,0,FSParams(si),fileList.Instrument,FS(si).Noise,FS(si).Base,NOISE); 
                    % Estimate SNR of selected peaks
                    PL(si).Noise = interp1(FS(si).f,FS(si).Noise,PL(si).f);
                    PL(si).Base = interp1(FS(si).f,FS(si).Base,PL(si).f);
                    PL(si).snr = (PL(si).h - PL(si).Base)./PL(si).Noise;
                    % flag peaks above SNR threshold
                    PL(si).flag = PL(si).snr > NOISE.MIN_SNR;
                end
            else
                % only pick peaks which have a height > dNthresh
                % (faster than picking all peaks then filtering)
                if NOISE.ESTIMATE
                    [PL(si).y, PL(si).h, PL(si).f, PL(si).res] = PeakDetectArea_0_3(FS(si).data,FS(si).f,peakMethod,FSParams(si).dNthresh,FSParams(si),fileList.Instrument,[],[],NOISE);
                else
                    [PL(si).y, PL(si).h, PL(si).f, PL(si).res] = PeakDetectArea_0_3(FS(si).data,FS(si).f,peakMethod,1,FSParams(si),fileList.Instrument,FS(si).Noise,FS(si).Base,NOISE);
                end
                % flag all picked peaks
                PL(si).flag = true(size(PL(si).y));
            end
        end
        % save peaks
        if DEBUG.SAVELOADPEAKS
            fprintf('Saving peaks to file %s\n',fName);
            save(fName,'PL');
        end
    end

    %% f. Calculate snr of each peak
    fprintf('Calculating peak SNR...\n');
    for si=1:numWindows
        % noise is peak height / standard deviation of data in signal-free
        % region
        if NOISE.ESTIMATE
            PL(si).snr = single(PL(si).h / FSParams(si).dNsd);
            PL(si).dNsd(1,1:length(PL(si).y)) = single(FSParams(si).dNsd); % Store with each peak the standard deviation of the noise in the SIM window
            PL(si).Base(1,1:length(PL(si).y)) = 0; % Set parameters from RAW noise info to zero since they are not used
            PL(si).Noise(1,1:length(PL(si).y)) = 0;
        else
            FSParams(si).dNthresh = 0; PL(si).dNsd(1,1:length(PL(si).y)) = 0; % Set parameters from manual noise estimation to zero. This way windows are stitched together at a place where there is a zero in the spectrum. See code around line 3315 
            if NOISE.INCLUDENOISEPEAKS_ON
                continue
            else
                PL(si).Noise = interp1(FS(si).f,FS(si).Noise,PL(si).f);
                PL(si).Base = interp1(FS(si).f,FS(si).Base,PL(si).f);
                PL(si).snr = single((PL(si).h - PL(si).Base) ./ PL(si).Noise);
            end
        end
    end
    % remove peak heights: not needed anymore
    PL = rmfield(PL,'h');
    
    %% g. Abundance correction = not used anymore. Step can be replaced with
    %normalization step to correct for intensity differences between SIM
    %windows.
    
    %% h. Remove "known" noise regions
    if NOISE.NOISEREMOVEALLKNOWN_ON
        nRegions = HIGHNOISEKNOWNREGIONS_MZ;
        fprintf('Removing peaks in noise regions...\n');
        for i=1:size(nRegions,1)
            fprintf('\t%.6f - %.6f m/z\n',nRegions(i,:));
        end
        fprintf('\n');
        for si=1:numWindows
            n = 0;    % number of data points removed from spectrum
            for j=1:size(nRegions,1)
                peaksmz = f2mz(PL(si).f,FSParams(si),fileList.Instrument); %Update for other machines
                idx = find(peaksmz < max(nRegions(j,:)) & peaksmz > min(nRegions(j,:)));
                % remove peak data
                sFields = fieldnames(PL);
                for i=1:length(sFields)
                    PL(si).(sFields{i})(idx) = [];
                end
                n = n+length(idx);
            end
            fprintf('Spectrum %d, removed %d data points\n',si,n);
        end
        fprintf('done\n');
    end
    
    %% i. Calibration
    
    % 1. Get list of internal calibrants 
    if CAL.ON && fileList.calOn(fi)
        %define what metabolites to use as internal calibrants and make sure
        %they're not in the exact mass list which we use for assessing the mass
        %accuracy of the spectrum
        if isempty(CAL.FILE)
            fprintf('No internal calibrant file - using external calibration (caution: results in sub-optimal stitching)');
            pause(2);
            calibrants = [];
            % change switches as not able to calibrate
            fileList.calOn(fi) = 0;
        else
            fprintf('Loading calibrant list file: %s...',CAL.FILE);
            try
                calfile = fopen(CAL.FILE); %24/04/2015 (JE) Replaced textread since this function will not be supported by matlab anymore soon
                calinfo = textscan(calfile,'%s %f','delimiter','\t');
                calibrants.ID = calinfo{1,1};
                calibrants.mz = calinfo{1,2};
                fclose(calfile);
                clear calfile; clear calinfo;
            catch
                error('Cannot load calibrant file %s',CAL.FILE);
            end
            fprintf('done\n');
            %remove any commented masses
            idx = find(strncmp(calibrants.ID,'%',1));
            calibrants.ID(idx) = [];
            calibrants.mz(idx) = [];
            %remove any duplicate m/z values and sort by increasing m/z
            [calibrants.mz, idx] = unique(calibrants.mz);
            calibrants.ID = calibrants.ID(idx);
            %remove any outside range of spectra
            mzMax = max([FSParams.mzEnd]);
            mzMin = min([FSParams.mzStart]);
            idx = find(mzMin<=calibrants.mz & calibrants.mz<=mzMax);
            calibrants.mz = calibrants.mz(idx).';
            calibrants.ID = calibrants.ID(idx).';
        end
    else
        % external calibration
        fprintf('External calibration only\n');
        calibrants = [];
        fileList.calOn(fi) = 0;
    end
    
    % 2. Internal recalibration 
    Pcal = FSParams;   % create new parameter set for calibrated SIM windows
    
    if CAL.ON && fileList.calOn(fi)
        % output information for user from calibration stage
        calOut = [];
        calOut.idx = [];
        calOut.sw = [];
        calOut.premz = [];
        calOut.y = [];
        calOut.snr = [];
        calOut.mode = [];
        calOut.postmz = [];
        calOut(1) = [];
        
        fprintf('Internally recalibrating SIM windows...\n');
        for si = 1:numWindows
            
            % initialise output struct for this SIM window
            sFields = fieldnames(calOut);
            calOutS = [];
            for i=1:length(sFields)
                calOutS.(sFields{i}) = [];
            end
            calOutS(1) = [];
            
            % potential calibrants
            idx = find(calibrants.mz >= Pcal(si).mzStart & calibrants.mz <= Pcal(si).mzEnd);
            nCals = length(idx);
            
            % are there any potential calibrants?
            if nCals > 0
                % yes: store in summary output
                for i=1:length(idx)
                    calOutS(i).idx = idx(i); % calibrant id
                    calOutS(i).sw = si;     % SIM window
                end
                potCals = calibrants.mz(idx);   % list potential calibrants (as stored in calOutS)
                
                % spectrum of interest: where flag set (detected peaks)
                idxF = find(PL(si).flag);
                
                % find closest peaks in spectrum
                mz = f2mz(PL(si).f(idxF),Pcal(si),fileList.Instrument); %Update for other machines
                intens = PL(si).y(idxF);
                [idxS,idxCal] = FindLargest_1_0(mz,potCals,CAL.MAX_PK_D,'ppm',intens);
                idxS = idxF(idxS);
                
                % store in output
                for i=1:length(idxCal)
                    mz = f2mz(PL(si).f(idxS(i)),Pcal(si),fileList.Instrument); %Update for other machines
                    calOutS(idxCal(i)).premz = mz;
                    calOutS(idxCal(i)).y = PL(si).y(idxS(i));
                    calOutS(idxCal(i)).snr = PL(si).snr(idxS(i));
                end
                
                % filter for peaks close enough to calibrants
                mz = f2mz(PL(si).f(idxS),Pcal(si),fileList.Instrument); %Update for other machines
                dmzppm = abs(mz-potCals(idxCal))./potCals(idxCal)*1e6;
                % update the indices for matching peaks
                idxS = idxS(dmzppm <= CAL.MAX_PK_D);
                idxCal = idxCal(dmzppm <= CAL.MAX_PK_D);
                
                % select those to use based on intensity
                idx = PL(si).snr(idxS) > CAL.MIN_SNR;
                % update the indices for matching peaks
                idxS = idxS(idx);
                idxCal = idxCal(idx);

                % update number of calibrants
                nCals = length(idxCal);
                
                % calibrants m/z
                calPoints = potCals(idxCal);
                % corresponding spectrum peaks
                PLcal = [];
                PLcal.y = PL(si).y(idxS);
                PLcal.f = PL(si).f(idxS);
            end
            fprintf('SIM window %2d - %2d calibrants found\n',si,nCals);
            
            % calibrated if possible
            if nCals == 0
                Pcal(si).icald = 0;    % No suitable calibrants found for this window
            else
                Pcal(si).icald = 1;
                %select actual calibration mode depending on how many calibration points there are
                switch CAL.MODE
                    case 'ab'
                        if (max(calPoints)-min(calPoints)) < (CAL.MIN_RANGE/100*(Pcal(si).mzEnd-Pcal(si).mzStart))
                            calMode = 'b';
                        elseif nCals >= 2
                            warning(['Using multi-parameter calibration over ',num2str((max(calPoints)-min(calPoints))/(Pcal(si).mzEnd-Pcal(si).mzStart)*100),'% but CAL.MIN_RANGE not determined accurately']);
                            calMode = 'ab';
                        else
                            calMode = 'b';
                        end
                    case 'abc'
                        switch(fileList.Instrument)
                            case{'ltqft'}
                                if (max(calPoints)-min(calPoints)) < (CAL.MIN_RANGE/100*(Pcal(si).mzEnd-Pcal(si).mzStart))
                                    calMode = 'b';
                                elseif nCals>=3
                                    warning(['Using multi-parameter calibration over ',num2str((max(calPoints)-min(calPoints))/(Pcal(si).mzEnd-Pcal(si).mzStart)*100),'% but CAL.MIN_RANGE not determined accurately']);
                                    calMode = 'abc';
                                elseif nCals>=2
                                    warning(['Using multi-parameter calibration over ',num2str((max(calPoints)-min(calPoints))/(Pcal(si).mzEnd-Pcal(si).mzStart)*100),'% but CAL.MIN_RANGE not determined accurately']);
                                    calMode = 'ab';
                                else
                                    calMode = 'b';
                                end
                            case{'qexactive','orbitrap'}
                                warning(['Three parameter calibration not available switching to two parameter mode']);
                                if (max(calPoints)-min(calPoints)) < (CAL.MIN_RANGE/100*(Pcal(si).mzEnd-Pcal(si).mzStart))
                                    calMode = 'b';
                                elseif nCals >= 2
                                    warning(['Using multi-parameter calibration over ',num2str((max(calPoints)-min(calPoints))/(Pcal(si).mzEnd-Pcal(si).mzStart)*100),'% but CAL.MIN_RANGE not determined accurately']);
                                    calMode = 'ab';
                                else
                                    calMode = 'b';
                                end
                        end
                end
                
                % get new calibration parameters
                [Pcal(si)] = Calibrate_1_0(PLcal, Pcal(si), calPoints, calMode, CAL.WEIGHTED_ON, fileList.Instrument); %Update for other machines
                
                % store information about the calibrants used
                mz = f2mz(PLcal.f,Pcal(si),fileList.Instrument); %Update for other machines
                for i=1:nCals
                    calOutS(idxCal(i)).mode = upper(calMode);
                    calOutS(idxCal(i)).postmz = mz(i);
                end
                
            end
            % update output
            calOut = [calOut calOutS];
        end
    else
        % no internal calibration
        for si=1:numWindows
            Pcal(si).icald = 0;
        end
    end
    
    % 3. External calibration
    % Use external calibration for windows not overlapping any internally
    % recalibrated windows - this is to reduce proliferation of errors through
    % uncalibrated windows
    % update calibration parameters Pcal
    
    intC = [Pcal.icald]; % window internally calibrated?
    
    for si=1:numWindows
        if ~intC(si)
            % not calibrated: how many overlapping int cal'd SIM windows?
            n = 0;
            for i=1:numWindows
                if intC(i) && (Pcal(i).mzStart < Pcal(si).mzEnd) && (Pcal(i).mzEnd > Pcal(si).mzStart)
                    n = n+1;
                end
            end
            if n == 0
                % use external calibration if overlapping no internally
                % recalibrated SIM windows
                Pcal(si).ecald = 1;
            else
                % don't calibrate if overlapping: align in next stage
                Pcal(si).ecald = 0;
            end
        else
            % internally calibrated: no need for external calibration
            Pcal(si).ecald = 0;
        end
    end
    
    % display summary of calibration % Change for new machines!
    for si=1:numWindows
        switch(fileList.Instrument)
            case{'ltqft'}
                fprintf('Window %2d, A=%.2f (d %+g), B=%10.5e (d %10.5e)\n',...
                    si, Pcal(si).A,Pcal(si).A-FSParams(si).A, ...
                    Pcal(si).B,Pcal(si).B-FSParams(si).B);
            case{'qexactive','orbitrap'}
                fprintf('Window %2d, B=%.2f (d %+g), C=%10.5e (d %10.5e)\n',...
                    si, Pcal(si).B,Pcal(si).B-FSParams(si).B, ...
                    Pcal(si).C,Pcal(si).C-FSParams(si).C);
            case{solarix}
        end
    end
    
    % j. Calibration via alignment with neighboring SIM windows.
    % align with adjacent internally or externally calibrated SIM windows
    
    % create new parameter set for aligned SIM windows
    Pali = Pcal;
    
    % output information for user from this stage: even if not aligning as
    % contains summary of calibration too
    aliOut = [];
    aliOut.sw = [];
    aliOut.cMode  = [];
    aliOut.npL = [];
    aliOut.npR = [];
    aliOut.mode = [];
    aliOut.preL = [];
    aliOut.preR = [];
    aliOut.dmzLoM = [];
    aliOut.dmzHi = [];
    aliOut.postL = [];
    aliOut.postR = [];
    aliOut(1) = [];
    
    if ALIGN.MZ_ALIGN_ON
        
        for si = 1:numWindows    %spectrum index
            
            % user output
            aliOut(si).sw = si;
            if Pcal(si).icald
                aliOut(si).cMode = 'i';
            elseif Pcal(si).ecald
                aliOut(si).cMode = 'e';
            else
                aliOut(si).cMode = 'a';
            end
            
            % calibrated?
            if Pcal(si).ecald || Pcal(si).icald
                % yes: no alignment
                Pali(si).alidL = 0; % aligned on left
                Pali(si).alidR = 0; % aligned on right
                Pali(si).m = [];
                Pali(si).c = [];
            else
                % no: align
                % overlapping SIM window indices
                % reference left & right (riL & riR)
                % "subject" spectrum being aligned (si)
                switch si
                    case 1
                        riL = 0;
                        riR = si+1;
                    case numWindows
                        riL = si-1;
                        riR = 0;
                    otherwise
                        riL = si-1;
                        riR = si+1;
                end
                fprintf('Aligning range %d-%d with\n',Pcal(si).mzStart,Pcal(si).mzEnd);
                if riL, fprintf('%d-%d on left\n',Pcal(riL).mzStart,Pcal(riL).mzEnd); end
                if riR, fprintf('%d-%d on right\n',Pcal(riR).mzStart,Pcal(riR).mzEnd); end
                
                % filter peaks in subject spectrum with SNR > ALIGN.MIN_SNR
                % AND flag set
                idx = find((PL(si).snr > ALIGN.MIN_SNR) & PL(si).flag);
                if riL, iSL = idx; end
                if riR, iSR = idx; end
                
                % filter peaks in reference spectrum with SNR > ALIGN.MIN_SNR
                
                if riL
                    iL = find(PL(riL).snr > ALIGN.MIN_SNR);
                end
                if riR
                    iR = find(PL(riR).snr > ALIGN.MIN_SNR);
                end
                
                
                % find common peaks and store
                if riL
                    % find closest with S
                    % using the original (uncalibrated) A and B for each
                    % spectrum gives the best indication of the closest peaks
                    mzSL = f2mz(PL(si).f(iSL),FSParams(si),fileList.Instrument);
                    mzL = f2mz(PL(riL).f(iL),FSParams(riL),fileList.Instrument);
                    [idxl, idxs] = FindClosest_0_3(mzL,mzSL,ALIGN.MAX_PK_D,'ppm');
                    if isempty(idxl)
                        % can't align with the left side!
                        riL = 0;
                    else
                        % back to calibrated spectra
                        mzSL = f2mz(PL(si).f(iSL),Pcal(si),fileList.Instrument);
                        mzL = f2mz(PL(riL).f(iL),Pcal(riL),fileList.Instrument);
                        % store common peaks to optimise alignment with
                        rL.mz = mzL(idxl);
                        rL.y = PL(riL).y(iL(idxl));
                        sL.mz = mzSL(idxs);
                        sL.y = PL(si).y(iSL(idxs));
                        sL.f = PL(si).f(iSL(idxs));
                        % store ref frequency with parameters homogenised to si
                        % so we can use the same parameters (Pcal(si)) for ref and si
                        rL.f = mz2f(rL.mz,Pcal(si),fileList.Instrument);
                    end
                    clear idxl idxs mzL mzSL
                end
                if riR
                    % find closest with S
                    % using the original (uncalibrated) A and B for each
                    % spectrum gives the best indication of the closest peaks
                    mzSR = f2mz(PL(si).f(iSR),FSParams(si),fileList.Instrument);
                    mzR = f2mz(PL(riR).f(iR),FSParams(riR),fileList.Instrument);
                    [idxr, idxs] = FindClosest_0_3(mzR,mzSR,ALIGN.MAX_PK_D,'ppm');
                    if isempty(idxr)
                        % can't align with the right side!
                        riR = 0;
                    else
                        % back to calibrated spectra
                        mzSR = f2mz(PL(si).f(iSR),Pcal(si),fileList.Instrument);
                        mzR = f2mz(PL(riR).f(iR),Pcal(riR),fileList.Instrument);
                        % store common peaks to optimise alignment with
                        rR.mz = mzR(idxr);
                        rR.y = PL(riR).y(iR(idxr));
                        sR.mz = mzSR(idxs);
                        sR.y = PL(si).y(iSR(idxs));
                        sR.f = PL(si).f(iSR(idxs));
                        % store ref frequency with parameters homogenised to si
                        % so we can use the same parameters (Pcal(si)) for ref and si
                        rR.f = mz2f(rR.mz,Pcal(si),fileList.Instrument);
                    end
                    clear idxr idxs mzR mzSR
                end
                
                % count number peaks on each side
                nL = 0;
                nR = 0;
                if riL, nL = length(rL.mz); end
                if riR, nR = length(rR.mz); end
                
                % user output
                aliOut(si).npL = nL;
                aliOut(si).npR = nR;
                
                % parameter output
                Pali(si).alidL = (nL>0);
                Pali(si).alidR = (nR>0);
                
                if (nL==0) && (nR==0)
                    %no alignment possible
                    fprintf('No peaks found in overlap to align, skipping...');
                else
                    % user output: RMS ppm error pre-alignment
                    if riL
                        % left
                        dmz = (rL.mz - sL.mz)./rL.mz *1e6; 
                        n = length(dmz);
                        aliOut(si).preL = norm(dmz)/sqrt(n);
                    end
                    if riR
                        % right
                        dmz = (rR.mz - sR.mz)./rR.mz *1e6;
                        n = length(dmz);
                        aliOut(si).preR = norm(dmz)/sqrt(n);
                    end
                    
                    % only use m scaling if there are multiple
                    % points and they are separated by at least
                    % ALIGN_M_MIN_RANGE(%) of the subject SIM window width
                    m = 1;    % initial m/z value scaling factor
                    c = 0;    % init  ial linear offset
                    if riL && riR
                        dmz = max([sR.mz sL.mz]) - min([sR.mz sL.mz]);
                    elseif riL
                        dmz = max(sL.mz) - min(sL.mz);
                    elseif riR
                        dmz = max(sR.mz) - min(sR.mz);
                    end
                    if ((nL+nR) < 2) || (dmz < ALIGN.M_MIN_RANGE/100*(Pcal(si).mzEnd-Pcal(si).mzStart))
                        fprintf('Using C parameter linear offset only\n');
                        params = c;
                        % update output summary file
                        aliOut(si).mode = 'o'; % zero order frequency: offset
                    else
                        % m = m/1e6 to increase sensitivity during fminsearch
                        params = [c;m*1e6]; %;a];
                        % update output summary file
                        aliOut(si).mode = 'so'; % zero and first order frequency: scale & offset
                    end
                    
                    % set spectra not being used in alignment to empty
                    if ~riR && riL
                        rR = [];
                        sR = [];
                    elseif ~riL && riR
                        rL = [];
                        sL = [];
                    end
                    
                    % optimise
                    paramTol=0.00000001;
                    covarTol = 1e99;
                    [params,fval] = fminsearch(@(params)AlignMinFunc_0_1(rL,sL,rR,sR,params), params,...
                        optimset('TolFun',covarTol,'TolX',paramTol)); %,'Display','iter'));
                    c = params(1);
                    if length(params) == 2
                        m = params(2) / 1e6;
                    else
                        m = 1;
                        %a = params(3);
                    end
                    
                    % new freq and m/z values
                    if riL
                        sL.f = m*sL.f+c;
                        sL.mz = f2mz(sL.f,Pcal(si),fileList.Instrument);
                    end
                    if riR
                        sR.f = m*sR.f+c;
                        sR.mz = f2mz(sR.f,Pcal(si),fileList.Instrument);
                    end
                    
                    % user output: RMS error post-alignment
                    if riL
                        % left
                        dmz = (rL.mz - sL.mz)./rL.mz *1e6;
                        n = length(dmz);
                        aliOut(si).postL = norm(dmz)/sqrt(n);
                    end
                    if riR
                        % right
                        dmz = (rR.mz - sR.mz)./rR.mz *1e6;
                        n = length(dmz);
                        aliOut(si).postR = norm(dmz)/sqrt(n);
                    end
                    
                    % store alignment parameters
                    Pali(si).m = m;
                    Pali(si).c = c;
                    %             Pali.a(si) = a;
                    
                end
                % clear temp spectra
                clear rR rL sR sL
            end
        end
        
        % apply adjustment parameters
        % fa = m*f + c
        for si=1:numWindows
            if Pali(si).alidL || Pali(si).alidR
                
                % display
                fprintf('Window %2d, m=%10.8e, c=%.2f\n',...
                    si, Pali(si).m,Pali(si).c);
                
                % apply parameters
                temp = PL(si).f;
                PL(si).f = Pali(si).m * PL(si).f + Pali(si).c;
                FS(si).f = Pali(si).m * FS(si).f + Pali(si).c;
                
                % user output: maximum potential change in m/z (occurs at high mass)
                dmz = f2mz(temp,Pali(si),fileList.Instrument) - f2mz(PL(si).f,Pali(si),fileList.Instrument);
                dmz = dmz ./ f2mz(temp,Pali(si),fileList.Instrument) * 1e6; % ppm
                [mz,idx] = max(abs(dmz));
                aliOut(si).dmzHi = dmz(idx);
                
            end
        end
        
    else
        fprintf('Skipping alignment\n');
        for si=1:numWindows
            % update parameters
            Pali(si).alidL = 0; % aligned on left
            Pali(si).alidR = 0; % aligned on right
            % user output
            aliOut(si).sw = si;
            aliOut(si).cMode = '';
            if Pcal(si).icald, aliOut(si).cMode = [aliOut(si).cMode, 'i']; end
            if Pcal(si).ecald, aliOut(si).cMode = [aliOut(si).cMode, 'e']; end
        end
    end
    
    % ####################################################################
    
    %% j. Homogenise all spectra to calibration parameters from SIM window 1
    
    Phom = Pali;    %create new parameter set for homogenized SIM windows
    
    fprintf('Homogenising calibration parameters');
    for si=1:numWindows
        switch(fileList.Instrument)
            case{'ltqft'}
                Phom(si).A = Pali(1).A;
                Phom(si).B = Pali(1).B;
                % 29/1/07 - final stitching does not need to use 3-term calibration, since it will add no value
                % at this point, and the only reason for stitching in the frequency
                % domain is to get the best peak shape for interpolation.
                Phom(si).C = 0;
            case{'qexactive','orbitrap'}
                Phom(si).A = 0;
                Phom(si).B = Pali(1).B;
                Phom(si).C = Pali(1).C;   
        end
        % data
        mz = f2mz(FS(si).f,Pali(si),fileList.Instrument);
        FS(si).f = mz2f(mz,Phom(si),fileList.Instrument);
        % peaks
        mz = f2mz(PL(si).f,Pali(si),fileList.Instrument);
        PL(si).f = mz2f(mz,Phom(si),fileList.Instrument);
        fprintf('.');
    end
    fprintf('done\n');
    
    %%  k. final stitch
    
    % new parameter set for final stitched spectrum
    % INCLUDE switch command for different machines. 
    Pfin.A = Phom(1).A;
    Pfin.B = Phom(1).B;
    Pfin.C = Phom(1).C;
    Pfin.T = Phom(1).T;  %since T=2^n/df and df and n are the same as for spec(1). % WHAT IS THIS?!
    Pfin.zFills = Phom(1).zFills;
    Pfin.mzStart = min([Phom.mzStart]);
    Pfin.mzEnd = max([Phom.mzEnd]);

    fprintf('Final stitch...\n');
    %currently spectra are interpolated over f points of FIRST spectrum!
    [FS, PL] = SegStitch_1_1(FS, PL, Phom, STITCH.NULL_REGION, fileList.Instrument);
    
    % sort peaks by ascending m/z
    peaksmz = f2mz(PL.f,Pfin,fileList.Instrument);
    [~, idx] = sort(peaksmz);
    sFields = fieldnames(PL);
    for i=1:length(sFields)
        PL.(sFields{i}) = PL.(sFields{i})(idx);
    end
    
    %% l. save results
    
    %save stitched spectrum as MATLAB file
    fprintf('Saving stitched spectrum...');
    if isempty(dir(SFSDir))
        fprintf('(making stitch dir)...');
        mkdir(SFSDir);
    end
    fileName = [SFSDir,opfname.stitchm];
    % recreate old data structures for compatibility
    SFS.data = FS.data;
    SFS.f = FS.f;
    SFS.peaksd = PL.y;
    SFS.peaksf = PL.f;
    SFS.peaksRes = PL.res;
    SFS.pSNR = PL.snr;
    SFS.pNoise = PL.Noise;
    SFS.pBase = PL.Base;
    SFS.dNsd = PL.dNsd;
    SFS.peaksFlag = PL.flag;
    SFSParams = Pfin;
    SFSParams.Instrument = fileList.Instrument; % Save instrument name for replicate filtering
    % convert some data to single to reduce storage
    SFS.data = single(SFS.data);
    SFS.peaksd = single(SFS.peaksd);
    SFS.peaksRes = single(SFS.peaksRes);
    SFS.pSNR = single(SFS.pSNR);
    % add extra info to parameters
    SFSParams.fileVersion = mfilename;
    % save
    save(fileName,'SFSParams','SFS');
    fprintf('done\n');
    clear SFS
    
    % save peak list as text file
    fprintf('Saving peak list text file...');
    if isempty(dir(SPLDir))
        fprintf('(making peaks dir)...');
        mkdir(SPLDir);
    end
    fileName = [SPLDir,'SPL_',fileList.Samples(1,fi).ID,'.txt'];
    % save file
    fid = fopen(fileName,'wt');
    if fid==-1
        error('Cannot create file %s',fileName);
    end
    peaksmz = f2mz(PL.f, Pfin,fileList.Instrument);
    if NOISE.NOISEFILT_ON
        % noise has been measured: output
        fprintf(fid,'M/Z\tINTENSITY\tSNR\tNON-NOISE_FLAG\n');
        for i=1:length(peaksmz)
           fprintf(fid,'%.7f\t%.6e\t%.2f\t%d\n',peaksmz(i),PL.y(i),PL.snr(i),PL.flag(i));
        end
    else
        % no noise measured
        fprintf(fid,'M/Z\tINTENSITY\n');
        for i=1:length(peaksmz)
            fprintf(fid,'%.7f\t%.6e\n',peaksmz(i),PL.y(i));
        end
    end
    fclose(fid);
    fprintf('done\n'); 
    
    % save internal calibration output
    if CAL.ON && fileList.calOn(fi)
        fprintf('Saving output from calibration...');
        fileName = [SFSDir,'calOutput','.txt'];
        if isempty(dir(fileName))
            % if file doesn't already exist, create and write headers
            fid = fopen(fileName,'w');
            fprintf(fid,['SAMPLE\tCAL\tEXACT_MASS\tSIM_WINDOW\tRANGE\tNEAREST_PK_M/Z\tPPM_ERROR\tINTENSITY\t',...
                'SNR\tCAL_MODE\tPOST_CAL_M/Z\tERROR_PPM\r\n']);
        else
            % if file exists, append
            fid = fopen(fileName,'a');
            fprintf(fid,['SAMPLE\tCAL\tEXACT_MASS\tSIM_WINDOW\tRANGE\tNEAREST_PK_M/Z\tPPM_ERROR\tINTENSITY\t',...
                'SNR\tCAL_MODE\tPOST_CAL_M/Z\tERROR_PPM\r\n']);
        end
        for i=1:length(calOut)
            % pre
            fprintf(fid,'%s\t%s\t%.8f\t%d\t%.2f-%.2f\t%.6f\t%.2f\t%.3e\t',...
                fileList.Samples(1,fi).ID,...
                calibrants.ID{calOut(i).idx},...
                calibrants.mz(calOut(i).idx),...
                calOut(i).sw,...
                FSParams(calOut(i).sw).mzStart,...
                FSParams(calOut(i).sw).mzEnd,...
                calOut(i).premz,...
                (calOut(i).premz-calibrants.mz(calOut(i).idx))./calibrants.mz(calOut(i).idx)*1e6,...
                calOut(i).y);
            % noise
            if NOISE.NOISEFILT_ON
                fprintf(fid,'%.2f\t',...
                    calOut(i).snr);
            else
                fprintf(fid,'u/k\t');
            end
            % post
            fprintf(fid,'%s\t%.6f\t%.2f\r\n',...
                calOut(i).mode,...
                calOut(i).postmz,...
                (calOut(i).postmz-calibrants.mz(calOut(i).idx))./calibrants.mz(calOut(i).idx)*1e6);
        end
        fclose(fid);
        fprintf('done\n');
    end
    
    % save alignment output
    fprintf('Saving output from m/z alignment...');
    fileName = align_file;
    if isempty(dir(fileName))
        % if file doesn't already exist, create and write headers
        fid = fopen(fileName,'wt');
        fprintf(fid,['SAMPLE\tSIM_WINDOW\tCAL_MODE\tNUM_PEAKS(L)\tNUM_PEAKS(R)\t',...
            'RMS_PPM_ERROR_PRE_(L)\tRMS_PPM_ERROR_PRE_(R)\tALI_MODE\tDMZ_PPM_HI\t',...
            'RMS_PPM_ERROR_POST_(L)\tRMS_PPM_ERROR_POST_(R)\n']);
    else
        % if file exists, append
        fid = fopen(fileName,'a');
        fprintf(fid,['SAMPLE\tSIM_WINDOW\tCAL_MODE\tNUM_PEAKS(L)\tNUM_PEAKS(R)\t',...
            'RMS_PPM_ERROR_PRE_(L)\tRMS_PPM_ERROR_PRE_(R)\tALI_MODE\tDMZ_PPM_HI\t',...
            'RMS_PPM_ERROR_POST_(L)\tRMS_PPM_ERROR_POST_(R)\n']);
    end
    for i=1:length(aliOut)
        fprintf(fid,'%s\t%2d\t%2s\t%3d\t%3d\t%.2f\t%.2f\t%.2s\t%.2f\t%.2f\t%.2f\n',...
            fileList.Samples(1,fi).ID,...
            aliOut(i).sw,...
            aliOut(i).cMode,...
            aliOut(i).npL,...
            aliOut(i).npR,...
            aliOut(i).preL,...
            aliOut(i).preR,...
            aliOut(i).mode,...
            aliOut(i).dmzHi,...
            aliOut(i).postL,...
            aliOut(i).postR);
    end
    fclose(fid);
    fprintf('done\n');
    
    if NOISE.ESTIMATE
        % save GetNoiseLevel output
        fprintf('Saving output from GetNoiseLevel...');
        fileName = noise_file;
        if isempty(dir(fileName))
            % if file doesn't already exist, create and write headers
            fid = fopen(fileName,'wt');
            fprintf(fid,['SAMPLE\tSIM_WINDOW\tNUMPOINTS\tSIGMA\tDNTHRESH\tMZSTART\tMZEND\n']);
        else
            % if file exists, append
            fid = fopen(fileName,'a');
            fprintf(fid,['SAMPLE\tSIM_WINDOW\tNUMPOINTS\tSIGMA\tDNTHRESH\tMZSTART\tMZEND\n']);
        end
        for i = 1:numWindows
            fprintf(fid,'%s\t%i\t%2d\t%2d\t%2d\t%3d\t%3d\n',...
                fileList.Samples(1,fi).ID,...
                i,...
                FSParams(i).dNnumpoints,...
                FSParams(i).dNsd,...
                FSParams(i).dNthresh,...
                FSParams(i).dNrange(1),...
                FSParams(i).dNrange(2));
        end
        fclose(fid);
        fprintf('done\n');
    end
end

% End of main loop

% Save all 

%create html output for Galaxy
%write Stitch html for Galaxy
fid = fopen(html_outfile_stitch, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>Stitched Frequency Spectra</b></h2></div>');
fprintf(fid,'<hr></hr>');

for i=1:length(stitchfile_html)
  	html_code = ['<div style="padding: 3px"><b><a href="', stitchfile_html{i},'">',stitchfile_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'<p>');
fprintf(fid,'');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);

%write peaklist html for Galaxy
fid = fopen(html_outfile_peaks, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>Stitched Peak Lists</b></h2></div>');
fprintf(fid,'<hr></hr>');
for i=1:length(peaklist_html)
  	html_code = ['<div style="padding: 3px"><b><a href="', peaklist_html{i},'">',peaklist_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'<p>');
fprintf(fid,'');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);
end


