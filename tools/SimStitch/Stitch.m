function Stitch_Galaxy(fileList, html_infile, html_indir, ListofCalibrants, csv_str, html_outfile_stitch, html_outdir_stitch, html_outfile_peaks, html_outdir_peaks, noise_file, align_file)
%
%   Stitch_Galaxy(config_file)
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
%           23,    sOptions.stitch.null_region.mzmin
%           24,    sOptions.stitch.null_region.mzmax
%           25,    sOptions.stitch.null_region.start
%           26,    sOptions.stitch.null_region.end
%           27,    sOptions.stitch.null_region.bound
%           28,    sOptions.segrange
%	    	29,    sOptions.calFile - This is now ListofCalibrants (see function input)
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


tmpName = tempname;
fid = fopen(tmpName,'w');
fprintf(fid, csv_str)
fclose(fid);

try 
    configs = csvread(tmpName);
	configs
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

%%%% Validating and assigning the full list of possible arguments with
%%%% reference to config_file -> configs


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
    %measure and filter noise? (off for reduced profile spectra)
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
    HIGHNOISEKNOWNREGIONS_MZ = [101.6 102.1; 101.32 101.42; 105.1 105.5; 74.05 74.2; 90.50 90.58; 116.37 116.5]; %DEFAULT
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
    CAL.WEIGHTED_ON = 1; %DEFAULT
    %is the calibration best fit weighted to the abundances?
end
if length(CAL.WEIGHTED_ON)>1   %VALIDATE 1
    sprintf('Error: config variable number 11. Too many values.')
    return
elseif isempty(CAL.WEIGHTED_ON)
    CAL.WEIGHTED_ON=0; %The line above ignores zeros, but was useful for checking extra values. 
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


temp = find(config_ind==23); %STITCH.NULL_REGION.MZMIN  ########################    23
if ~isempty(temp)
    STITCH.NULL_REGION.MZMIN = configs(temp,2:max(find(configs(temp,1:end))));
else    
    STITCH.NULL_REGION.MZMIN = 0; %DEFAULT
    %edge effect to remove from start of first SIM window
end
if length(STITCH.NULL_REGION.MZMIN)>1   %VALIDATE 1
    sprintf('Error: config variable number 23. Too many values.')
    return
elseif isempty(STITCH.NULL_REGION.MZMIN)
    STITCH.NULL_REGION.MZMIN=0; %The line above ignores zeros, but was useful for checking extra values. 
end
    
    
temp = find(config_ind==24); %STITCH.NULL_REGION.MZMAX  ########################    24
if ~isempty(temp)
    STITCH.NULL_REGION.MZMAX = configs(temp,2:max(find(configs(temp,1:end))));
else    
    STITCH.NULL_REGION.MZMAX = 0; %DEFAULT
     %edge effect to remove from end of last SIM window
end
if length(STITCH.NULL_REGION.MZMAX)>1   %VALIDATE 1
    sprintf('Error: config variable number 24. Too many values.')
    return
elseif isempty(STITCH.NULL_REGION.MZMAX)
    STITCH.NULL_REGION.MZMAX=0; %The line above ignores zeros, but was useful for checking extra values. 
end


temp = find(config_ind==25); %STITCH.NULL_REGION.START  ########################    25
if ~isempty(temp)
    STITCH.NULL_REGION.START = configs(temp,2:max(find(configs(temp,1:end))));
else    
    STITCH.NULL_REGION.START = 15; %DEFAULT
    %null or dead region (in m/z) from start of window
end
if length(STITCH.NULL_REGION.START)>1   %VALIDATE 1
    sprintf('Error: config variable number 25. Too many values.')
    return
elseif isempty(STITCH.NULL_REGION.START)
    STITCH.NULL_REGION.START=0; %The line above ignores zeros, but was useful for checking extra values. 
end


temp = find(config_ind==26); %STITCH.NULL_REGION.END    ########################    26
if ~isempty(temp)
    STITCH.NULL_REGION.END = configs(temp,2:max(find(configs(temp,1:end))));
else    
    STITCH.NULL_REGION.END = 15; %DEFAULT
    %null or dead region (in m/z) from end of window
end
if length(STITCH.NULL_REGION.END)>1   %VALIDATE 1
    sprintf('Error: config variable number 26. Too many values.')
    return
elseif isempty(STITCH.NULL_REGION.END)
    STITCH.NULL_REGION.END=0; %The line above ignores zeros, but was useful for checking extra values. 
end



temp = find(config_ind==27); %STITCH.NULL_REGION.BOUND      ####################    27
if ~isempty(temp)
    STITCH.NULL_REGION.BOUND = configs(temp,2:max(find(configs(temp,1:end))));
else    
    STITCH.NULL_REGION.BOUND = [70 2000]; %DEFAULT
    %                          boundary for each dead region (midpoint of OVERLAP of SIM windows)
    %                          It is possible to define different size edges for different regions of the spectrum:
    %                          [a b c d]: a < mz <= b
    %                                     b < mz <= c
    %                                     c < mz <= d;  mz is the midpoint of overlap
end

if isempty(STITCH.NULL_REGION.BOUND)
    STITCH.NULL_REGION.BOUND= [0 0] ; %The line above ignores zeros, but was useful for checking extra values. 
    % NB! THESE VALUES WOULD PROBABLY CRASH THE SYSTEM - OR AT LEAST MEAN
    % NO STITCHING IS DONE
end
    

% #####################ADDED FOR MSMS PROCESSING - WINDOW SELECT
temp = find(config_ind==28); %SEG.RANGE    ########################    28
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

% ################## calFile now part of sOptions
%temp = find(config_ind==29); %CAL.FILE      ####################    27
%if ~isempty(temp)
CAL.FILE = ListofCalibrants;
%else    
%    CAL.FILE = ''; %DEFAULT - no list of internal calibrants
%    
%end

% ---DEBUG--- %
DEBUG.QUICKPEAKDETECT = 0;    % speed up peak-peaking for debug only (not accurate)
DEBUG.SAVELOADPEAKS = 0;      % save or load intermediate list of detected peaks in _peaks.mat file?

%ex-global parameters
DISPLAY_HDRS = 0;   %display file header info
ZFILLS = 1;        %the number of zero-fills to be applied to transients
X_ZFILLS = 1;      %zero-fills assumed to have been applied by xCalibur for .raw file spectra

% ---FILE--- %
BLANKSTR = 'blank';                 % group name for blank samples


% ---NOISE FILTERING--- %

% ---INTERNAL CALIBRATION--- %

% ---ALIGNMENT--- %

% ---STITCHING--- %

STITCH.ON = 1;

% ---TIC VARS ---%
TICfilter = [];
avTICs= [];
RAWname = {};





% ########################################################################
%%  GET LIST OF FILES
% ########################################################################

%extract the fileList information from the xml file
fileList = import_filelist_xml(fileList);

%adjust fileList directories to match Galaxy temporary html file paths
fileList.specDir = [html_indir,'/']; %this converts the path given by Galaxy for the files output by ProcessTransients, adding a forward slash (may need changed for other OS)
fileList.stitchDir = [html_outdir_stitch,'/']; %this updates the fileList directory instructions for files to be output by Stitch (Galaxy output no 1)
fileList.peaksDir = [html_outdir_peaks, '/']; %this updates the fileList directory instructions for peak files to be output by Stitch (Galaxy output no 2)

%make output directories
dFold = fileList.stitchDir;
[sdir, mess, messid] = mkdir(dFold);
if ~sdir, error(['Cannot make stitched spectra directory ',dFold,': ',messid]); end

dFold = fileList.peaksDir;
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

% ########################################################################
%% MAIN LOOP 
% ########################################################################

stitchfile_html = {}; %these are place holders to catch the names of output files to be put into ...
peaklist_html = {};  %...html files for Galaxy to use. See collation at end of main loop. 

rawCount = length(fileList.rawFull);
% loop each raw file
for fi=1:rawCount
    
    % output filenames
    opfname.stitchm = ['specOut_',fileList.setUID,'_',fileList.rawShort{fi},'.mat'];
    
    %if we want to save the result, check if already stitched
    fprintf('\nFile: %s\n',fileList.rawShort{fi});
    fprintf('Looking for stitched file...');
    fileName = [fileList.stitchDir,opfname.stitchm];
    fileName = dir(fileName);
    if ~isempty(fileName), fprintf('stitched file exists, skipping\n'); continue;
    else fprintf('not found\n'); end
    
    % label whether or not this is a "blank" sample
    if strncmpi(fileList.spec(fi).sampleID,BLANKSTR,length(BLANKSTR))
        fileList.isBlank(fi) = true;
    else
        fileList.isBlank(fi) = false;
    end
    
    % label whether or not this sample should be internally recalibrated
    % (global cal is on and this isn't a blank or it is and we want to
    % calibrate blanks)
    if CAL.ON && (~fileList.isBlank(fi) || (fileList.isBlank(fi) && CAL.BLANK_ON))
        fileList.calOn(fi) = true;
    else
        fileList.calOn(fi) = false;
    end
    
    % ####################################################################
    %%  GET (UNCALIBRATED) FREQUENCY SPECTRA
    % ####################################################################
    
    if FILES.USERAWONLY_ON
        % use reduced profile data from .raw files only
        options.getSpec = 1;
        options.sumScans = 1;
        fprintf('Loading reduced profile data from %s.raw\n',fileList.rawShort{fi});
        % get spectral data and parameters from raw file
        [SD,P] = GetRawProfileFS_v2(fileList,options,fi,DISPLAY_HDRS,X_ZFILLS);
        NOISE.NOISEFILT_ON = 0;   % turn off noise filtering
    else
        fileName = [fileList.specDir,'spec',fileList.setUID,'_',fileList.rawShort{fi},'.mat'];
        fprintf('Loading pre-process spectral file: %s...',fileName);
        try
            load(fileName);
        catch
            fprintf('not found!\n');
            rethrow(lasterror);
        end
        % change field names
        P = specParams;
        SD = spec;
        for si=1:length(P)
            SD(si).y = spec(si).data;
        end
        SD = rmfield(SD,'data');
        clear spec specParams
        fprintf('done\n');
    end
    
    %the number of spectra (segments or windows)
    
    
    % ##########################ADDED FOR MSMS PROCESSING
    if ~isempty(SEG.RANGE) %reduce P and SD to just the chosen range
        if SEG.RANGE(2)>length(P) %check that range is within bounds.
            sprintf('Error: SEG RANGE (2) is greater than total number of windows!')
            return
        else
            P = P(SEG.RANGE(1):SEG.RANGE(2));
            SD = SD(SEG.RANGE(1):SEG.RANGE(2));
        end
    end
        
    
    numSpectra = length(P);
    
    %check data and freq's same length
    for si=1:numSpectra
        if length(SD(si).y)~=length(SD(si).f), error('Spectral intensity and f size mismatch'); end
    end
    
    % ####################################################################
    %% GET TIC
    % ####################################################################
    TIC = [P.TIC];
    
    
    % ####################################################################
    %%  CHECK FOR UNEXPECTED PARAMETERS
    % ####################################################################
    
    fprintf('Checking parameters...');
    % VT always 1
    if ~FILES.USERAWONLY_ON
        
        if ~isempty(find([P.VT]~=1, 1))
            error('Unexpected changing or non-unity VT value');
        end
    else
        %ignore VT parameter as it is not produced in the RAW file
        %parameters for some reason.
        [P.VT] = deal(1); % create VT and set to 1 for each spec.
    end
    
    % C always zero
    if ~isempty(find([P.C]~=0, 1))
        error('Unexpected non-zero C parameter in input files');
    end
    fprintf('done\n');
    
    % ####################################################################
    %%  SAMPLING FREQUENCY CORRECTION
    % ####################################################################
    % Correct the calculated spectra for SIM windows with a reduced sampling
    % frequency: see notebook P145.  Windows with lower sampling frequency
    % should have spectral y-values increased proportionally to correct for
    % reduced data points.
    
   
    %###################################################################
    %#################### AM REMOVING THIS FOR RAWONLY #################
    %###################################################################
    if ~FILES.USERAWONLY_ON
        fprintf('Correcting for varying sampling frequency...');
        % find maximum sampling frequency
        fs_max = max([P.BW]) * 2;
        
        % correct each SIM window
        for si=1:numSpectra
            % gamma is multiplier to achieve maximum fs
            gamma = fs_max / (P(si).BW * 2);
            SD(si).y = SD(si).y * gamma;
        end
        fprintf('done\n');
    end
    % ####################################################################
    %%  GET LIST OF INTERNAL CALIBRANTS
    % ####################################################################
    
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
                [calibrants.ID,calibrants.mz] = textread(CAL.FILE,'%s %f','delimiter','\t');
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
            mzMax = max([P.mzEnd]);
            mzMin = min([P.mzStart]);
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
    
    % ####################################################################
    %%  ESTIMATE NOISE IN SPECTRA
    % ####################################################################
    if NOISE.NOISEFILT_ON
        %first we determine the noise level of each spectrum by fitting an area
        %with estimated no signal to the Rayleigh distribution
        fprintf('Measuring noise\n');
        
        for si = 1:numSpectra
            % regions of spectrum (in Hz) to NOT use for noise measurement
            exRegion = mz2f(HIGHNOISEKNOWNREGIONS_MZ,P(si));
            % get data rms noise, sigma of Rayleigh fit, and freq range over which it is measured
            [nRMS,nSigma,nfRange, numPoints] = GetNoiseLevel_0_4(SD(si).y,SD(si).f,exRegion);
            P(si).dNrange = f2mz(fliplr(nfRange),P(si));
            % standard deviation of noise (in real, imaginary components of FT,
            % assumed equal and zero-mean Gaussian, estimated from Rayleigh
            % distribution fit of pure noise in signal-free region of
            % amplitude-mode spectrum).
            P(si).dNsd = nSigma;
            % peak height threshold
            P(si).dNthresh = NOISE.MIN_SNR * P(si).dNsd;
            fprintf('dNthresh=%.2f (%.2f x noise sd, range %.3f-%.2f)\n', ...
                P(si).dNthresh,NOISE.MIN_SNR,P(si).dNrange(1),P(si).dNrange(2));
            
            P(si).dNnumpoints = numPoints; %added by RLD 14/02/2013
        end
    else
        % no noise filtering (ie stitching reduced-profile spectra)
        for si = 1:numSpectra
            % set threshold to zero
            P(si).dNthresh = 0;
            P(si).dNsd = 0;
            P(si).dNnumpoints = 0;
            P(si).dNrange = [0 0];
            
        end
        % force picking of non-noise peaks only (which is all peaks since
        % threshold is set to 0, but then forces flag to be set on all peaks)
        NOISE.INCLUDENOISEPEAKS_ON = 0;
    end
    
    % ####################################################################
    %% THRESHOLD DATA
    % ####################################################################
    % only run if in "quick" mode where not all peaks are required: speeds up
    % peak-picking
    if NOISE.NOISEFILT_ON && ~NOISE.INCLUDENOISEPEAKS_ON
        % apply hard threshold to zero all data below noise
        % 01/Aug/07: don't threshold any data point adjacent to maxima ABOVE the threshold;
        % these will be used to interpolate to estimate the peak area
        for si = 1:numSpectra
            peaksn = maxima(SD(si).y,length(SD(si).y)) + 1;
            % only interested in peaks with maximum data point above the nf
            peaksn = peaksn(SD(si).y(peaksn) > P(si).dNthresh);
            fprintf('Window %d, number peaks > thresh = %d\n',si,length(peaksn));
            % and the points either side
            peaksn = unique([peaksn peaksn-1 peaksn+1]);
            % all points above nf
            idx = find((SD(si).y - P(si).dNthresh) >= 0);
            % above nf OR adjacent to peak
            idx = union(idx, peaksn);
            % all other points to zero
            idx = setdiff(1:length(SD(si).y), idx);
            SD(si).y(idx) = 0;
        end
    end
    
    % ####################################################################
    %%  DETECT, MEASURE AND FLAG PEAKS
    % ####################################################################
    % check for existing _peaks file (debug, saves time to dump peaks)
    fName = [fileList.peaksDir,fileList.setUID,fileList.rawShort{fi},'_peaks.mat'];
    files = dir(fName);
    if DEBUG.SAVELOADPEAKS && ~isempty(files)
        fprintf('Loading peaks from file %s\n',fName);
        load(fName);
    else
        if DEBUG.QUICKPEAKDETECT
            peakMethod = 'centroid';
        else
            peakMethod = 'KCEinterp';
        end
        fprintf('Finding peaks:\n');
        for si = 1:numSpectra
            fprintf('SIM window %d...',si);
            if NOISE.INCLUDENOISEPEAKS_ON
                % pick all peaks, including all noise peaks
                [SP(si).y, SP(si).h, SP(si).f, SP(si).res] = PeakDetectArea_0_3(SD(si).y,SD(si).f,peakMethod,0);
                % flag peaks with height > threshold
                SP(si).flag = SP(si).h > P(si).dNthresh;
            else
                % only pick peaks which have a height > dNthresh
                % (faster than picking all peaks then filtering)
                [SP(si).y, SP(si).h, SP(si).f, SP(si).res] = PeakDetectArea_0_3(SD(si).y,SD(si).f,peakMethod,P(si).dNthresh);
                % flag all picked peaks
                SP(si).flag = true(size(SP(si).y));
            end
        end
        % save peaks
        if DEBUG.SAVELOADPEAKS
            fprintf('Saving peaks to file %s\n',fName);
            save(fName,'SP');
        end
    end
    
    % ####################################################################
    %% CALCULATE SNR OF PEAKS
    % ####################################################################
    if NOISE.NOISEFILT_ON
        fprintf('Calculating peak SNR...\n');
        for si=1:numSpectra
            % noise is peak height / standard deviation of data in signal-free
            % region
            SP(si).snr = single(SP(si).h / P(si).dNsd);
        end
    else
        % can't calculate snr
        for si = 1:numSpectra
            SP(si).snr = single(zeros(size(SP(si).f)));
            %SP(si).snr = single(zeros(size(SP(si).f)));
        end
    end
    % remove peak heights: not needed anymore
    SP = rmfield(SP,'h');
    
    % ####################################################################
    %% STORE NOISE LOCAL TO EACH PEAK
    % ####################################################################
    % store with each peak the standard deviation of the noise in the SIM window
    for si=1:numSpectra
        SP(si).dNsd(1,1:length(SP(si).y)) = single(P(si).dNsd);
    end
    
    % ####################################################################
    %%  ABUNDANCE CORRECTION
    % ####################################################################
    % correct abundances to minimise residual difference between overlapping windows
    % parameterised correction is applied to each window before any further processing
    % correction is applied to PEAKS ONLY (ie spectral intensities will be UNcorrected)
    % Noise is assumed to be white (ie is not corrected)
    
    if STITCH.ON && ALIGN.INTENSITY_CORRECT_ON
        fprintf('Correcting intensity measurements in SIM windows...\n');
        fprintf('Sum of squares of log of flagged peak intensity error in overlap:\n');
        for si=1:numSpectra
            % reference window index
            if si>1, ri=si-1; else ri=0; end
            if ri
                % only use flagged peaks
                idx = SP(ri).flag;
                rd = SP(ri).y(idx);
                rmz = f2mz(SP(ri).f(idx),P(ri),rd);
                idx = SP(si).flag;
                sd = SP(si).y(idx);
                smz = f2mz(SP(si).f(idx),P(si),sd);
                [idx1 idx2] = FindClosest_0_3(rmz,smz,ALIGN.MAX_PK_D,'ppm');
                if ~isempty(idx1)
                    ss1 = sum((log10(rd(idx1)) - log10(sd(idx2))).^2);
                    fprintf('SS(log error), SIM window %3d, before=%.2f, ',si,ss1);
                end
            end
            % correction
            m = mean([P(si).mzEnd P(si).mzStart]);
            gamma = CAL.AC_M*m + CAL.AC_C;
            peaksmz = f2mz(SP(si).f,P(si),SP(si).y);
            d = peaksmz - P(si).mzStart;
            alpha = gamma*(d-CAL.AC_P) + 1;
            SP(si).y = SP(si).y ./ alpha;
            % sum of squares error in overlap after correction
            if ri
                % only use flagged peaks
                idx = SP(ri).flag;
                rd = SP(ri).y(idx);
                rmz = f2mz(SP(ri).f(idx),P(ri),rd);
                idx = SP(si).flag;
                sd = SP(si).y(idx);
                smz = f2mz(SP(si).f(idx),P(si),sd);
                [idx1 idx2] = FindClosest_0_3(rmz,smz,ALIGN.MAX_PK_D,'ppm');
                if ~isempty(idx1)
                    ss = sum((log10(rd(idx1)) - log10(sd(idx2))).^2);
                    fprintf('after=%.2f, d=%.2f\n',ss,ss-ss1);
                end
            end
        end
    end
    
    % ####################################################################
    %% REMOVE "KNOWN" REGIONS OF HIGH NOISE
    % ####################################################################
    
    if NOISE.NOISEREMOVEALLKNOWN_ON
        nRegions = HIGHNOISEKNOWNREGIONS_MZ;
        fprintf('Removing peaks in noise regions...\n');
        for i=1:size(nRegions,1)
            fprintf('\t%.6f - %.6f m/z\n',nRegions(i,:));
        end
        fprintf('\n');
        for si=1:numSpectra
            n = 0;    % number of data points removed from spectrum
            for j=1:size(nRegions,1)
                peaksmz = f2mz(SP(si).f,P(si));
                idx = find(peaksmz < max(nRegions(j,:)) & peaksmz > min(nRegions(j,:)));
                % remove peak data
                sFields = fieldnames(SP);
                for i=1:length(sFields)
                    SP(si).(sFields{i})(idx) = [];
                end
                n = n+length(idx);
            end
            fprintf('Spectrum %d, removed %d data points\n',si,n);
        end
        fprintf('done\n');
    end
    
    % ####################################################################
    %%  INTERNALLY RECALIBRATE
    % ####################################################################
    
    Pcal = P;   % create new parameter set for calibrated SIM windows
    
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
        for si = 1:numSpectra
            
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
                
                % spectrum of interest: where flag set
                idxF = find(SP(si).flag);
                
                % find closest peaks in spectrum
                mz = f2mz(SP(si).f(idxF),Pcal(si));
                intens = SP(si).y(idxF);
                [idxS,idxCal] = FindLargest_1_0(mz,potCals,CAL.MAX_PK_D,'ppm',intens);
                idxS = idxF(idxS);
                
                % store in output
                for i=1:length(idxCal)
                    mz = f2mz(SP(si).f(idxS(i)),Pcal(si));
                    calOutS(idxCal(i)).premz = mz;
                    calOutS(idxCal(i)).y = SP(si).y(idxS(i));
                    calOutS(idxCal(i)).snr = SP(si).snr(idxS(i));
                end
                
                % filter for peaks close enough to calibrants
                mz = f2mz(SP(si).f(idxS),Pcal(si));
                dmzppm = abs(mz-potCals(idxCal))./potCals(idxCal)*1e6;
                % update the indices for matching peaks
                idxS = idxS(dmzppm <= CAL.MAX_PK_D);
                idxCal = idxCal(dmzppm <= CAL.MAX_PK_D);
                
                % select those to use based on intensity
                idx = SP(si).snr(idxS) > CAL.MIN_SNR;
                % update the indices for matching peaks
                idxS = idxS(idx);
                idxCal = idxCal(idx);
                
                % update number of calibrants
                nCals = length(idxCal);
                
                % calibrants m/z
                calPoints = potCals(idxCal);
                % corresponding spectrum peaks
                SPcal = [];
                SPcal.y = SP(si).y(idxS);
                SPcal.f = SP(si).f(idxS);
            end
            fprintf('SIM window %2d - %2d calibrants found\n',si,nCals);
            
            % calibrated if possible
            if nCals == 0
                Pcal(si).icald = 0;    % internally recalibrated?
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
                end
                
                % get new calibration parameters
                [Pcal(si)] = Calibrate_1_0(SPcal, Pcal(si), calPoints, calMode, CAL.WEIGHTED_ON);
                
                % store information about the calibrants used
                mz = f2mz(SPcal.f,Pcal(si));
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
        for si=1:numSpectra
            Pcal(si).icald = 0;
        end
    end
    
    % ####################################################################
    %% EXTERNAL CALIBRATION
    % ####################################################################
    % Use external calibration for windows not overlapping any internally
    % recalibrated windows - this is to reduce proliferation of errors through
    % uncalibrated windows
    
    % update calibration parameters Pcal
    
    intC = [Pcal.icald]; % window internally calibrated?
    
    for si=1:numSpectra
        if ~intC(si)
            % not calibrated: how many overlapping int cal'd SIM windows?
            n = 0;
            for i=1:numSpectra
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
    
    % display summary of calibration
    for si=1:numSpectra
        fprintf('Window %2d, A=%.2f (d %+g), B=%10.5e (d %10.5e)\n',...
            si, Pcal(si).A,Pcal(si).A-P(si).A, ...
            Pcal(si).B,Pcal(si).B-P(si).B);
    end
    
    % ####################################################################
    %%  ALIGN SIM WINDOWS
    % ####################################################################
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
        
        for si = 1:numSpectra    %spectrum index
            
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
                    case numSpectra
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
                idx = find((SP(si).snr > ALIGN.MIN_SNR) & SP(si).flag);
                if riL, iSL = idx; end
                if riR, iSR = idx; end
                
                % filter peaks in reference spectrum with SNR > ALIGN.MIN_SNR
                if riL
                    iL = find(SP(riL).snr > ALIGN.MIN_SNR);
                end
                if riR
                    iR = find(SP(riR).snr > ALIGN.MIN_SNR);
                end
                
                % find common peaks and store
                if riL
                    % find closest with S
                    % using the original (uncalibrated) A and B for each
                    % spectrum gives the best indication of the closest peaks
                    mzSL = f2mz(SP(si).f(iSL),P(si));
                    mzL = f2mz(SP(riL).f(iL),P(riL));
                    [idxl idxs] = FindClosest_0_3(mzL,mzSL,ALIGN.MAX_PK_D,'ppm');
                    if isempty(idxl)
                        % can't align with the left side!
                        riL = 0;
                    else
                        % back to calibrated spectra
                        mzSL = f2mz(SP(si).f(iSL),Pcal(si));
                        mzL = f2mz(SP(riL).f(iL),Pcal(riL));
                        % store common peaks to optimise alignment with
                        rL.mz = mzL(idxl);
                        rL.y = SP(riL).y(iL(idxl));
                        sL.mz = mzSL(idxs);
                        sL.y = SP(si).y(iSL(idxs));
                        sL.f = SP(si).f(iSL(idxs));
                        % store ref frequency with parameters homogenised to si
                        % so we can use the same parameters (Pcal(si)) for ref and si
                        rL.f = mz2f(rL.mz,Pcal(si));
                    end
                    clear idxl idxs mzL mzSL
                end
                if riR
                    % find closest with S
                    % using the original (uncalibrated) A and B for each
                    % spectrum gives the best indication of the closest peaks
                    mzSR = f2mz(SP(si).f(iSR),P(si));
                    mzR = f2mz(SP(riR).f(iR),P(riR));
                    [idxr idxs] = FindClosest_0_3(mzR,mzSR,ALIGN.MAX_PK_D,'ppm');
                    if isempty(idxr)
                        % can't align with the right side!
                        riR = 0;
                    else
                        % back to calibrated spectra
                        mzSR = f2mz(SP(si).f(iSR),Pcal(si));
                        mzR = f2mz(SP(riR).f(iR),Pcal(riR));
                        % store common peaks to optimise alignment with
                        rR.mz = mzR(idxr);
                        rR.y = SP(riR).y(iR(idxr));
                        sR.mz = mzSR(idxs);
                        sR.y = SP(si).y(iSR(idxs));
                        sR.f = SP(si).f(iSR(idxs));
                        % store ref frequency with parameters homogenised to si
                        % so we can use the same parameters (Pcal(si)) for ref and si
                        rR.f = mz2f(rR.mz,Pcal(si));
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
                        dmz = (rL.mz - sL.mz)./sL.mz *1e6;
                        n = length(dmz);
                        aliOut(si).preL = norm(dmz)/sqrt(n);
                    end
                    if riR
                        % right
                        dmz = (rR.mz - sR.mz)./sR.mz *1e6;
                        n = length(dmz);
                        aliOut(si).preR = norm(dmz)/sqrt(n);
                    end
                    
                    % only use m scaling if there are multiple
                    % points and they are separated by at least
                    % ALIGN_M_MIN_RANGE(%) of the subject SIM window width
                    m = 1;    % initial m/z value scaling factor
                    c = 0;    % initial linear offset
                    if riL & riR
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
                    if ~riR & riL
                        rR = [];
                        sR = [];
                    elseif ~riL & riR
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
                        sL.mz = f2mz(sL.f,Pcal(si));
                    end
                    if riR
                        sR.f = m*sR.f+c;
                        sR.mz = f2mz(sR.f,Pcal(si));
                    end
                    
                    % user output: RMS error post-alignment
                    if riL
                        % left
                        dmz = (rL.mz - sL.mz)./sL.mz *1e6;
                        n = length(dmz);
                        aliOut(si).postL = norm(dmz)/sqrt(n);
                    end
                    if riR
                        % right
                        dmz = (rR.mz - sR.mz)./sR.mz *1e6;
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
        % fa = m*f + c + a*f.^2;
        for si=1:numSpectra
            if Pali(si).alidL || Pali(si).alidR
                
                % display
                fprintf('Window %2d, m=%10.8e, c=%.2f\n',...
                    si, Pali(si).m,Pali(si).c);
                
                % apply parameters
                temp = SP(si).f;
                SP(si).f = Pali(si).m * SP(si).f + Pali(si).c;
                SD(si).f = Pali(si).m * SD(si).f + Pali(si).c;
                
                % user output: maximum potential change in m/z (occurs at high mass)
                dmz = f2mz(temp,Pali(si)) - f2mz(SP(si).f,Pali(si));
                dmz = dmz ./ f2mz(temp,Pali(si)) * 1e6; % ppm
                [mz,idx] = max(abs(dmz));
                aliOut(si).dmzHi = dmz(idx);
                
            end
        end
        
    else
        fprintf('Skipping alignment\n');
        for si=1:numSpectra
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
    %%  HOMOGENISE ALL SPECTRA TO PARAMETERS FROM SIM WINDOW 1
    % ####################################################################
    
    Phom = Pali;    %create new parameter set for homogenized SIM windows
    
    fprintf('Homogenising parameters');
    for si=1:numSpectra
        Phom(si).A = Pali(1).A;
        Phom(si).B = Pali(1).B;
        % 29/1/07 - final stitching does not need to use 3-term calibration, since it will add no value
        % at this point, and the only reason for stitching in the frequency
        % domain is to get the best peak shape for interpolation.
        Phom(si).C = 0;
        % data
        mz = f2mz(SD(si).f,Pali(si));
        SD(si).f = mz2f(mz,Phom(si));
        % peaks
        mz = f2mz(SP(si).f,Pali(si));
        SP(si).f = mz2f(mz,Phom(si));
        fprintf('.');
    end
    fprintf('done\n');
    
    % ####################################################################
    %%  FINAL STITCH
    % ####################################################################
    
    % new parameter set for final stitched spectrum
    Pfin.A = Phom(1).A;
    Pfin.B = Phom(1).B;
    Pfin.C = Phom(1).C;
    Pfin.T = Phom(1).T;  %since T=2^n/df and df and n are the same as for spec(1).
    Pfin.zFills = Phom(1).zFills;
    Pfin.mzStart = min([Phom.mzStart]);
    Pfin.mzEnd = max([Phom.mzEnd]);
    
    fprintf('Final stitch...\n');
    %currently spectra are interpolated over f points of FIRST spectrum!
    [SD, SP] = SegStitch_1_1(SD, SP, Phom, STITCH.NULL_REGION);
    
    % sort peaks by ascending m/z
    peaksmz = f2mz(SP.f,Pfin);
    [peaksmz, idx] = sort(peaksmz);
    sFields = fieldnames(SP);
    for i=1:length(sFields)
        SP.(sFields{i}) = SP.(sFields{i})(idx);
    end
    
    % ####################################################################
    %%  SAVE RESULT
    % ####################################################################
    
    %save stitched spectrum as MATLAB file
    fprintf('Saving stitched spectrum...');
    if isempty(dir(fileList.stitchDir))
        fprintf('(making stitch dir)...');
        mkdir(fileList.stitchDir);
    end
    fileName = [fileList.stitchDir,opfname.stitchm];
    % recreate old data structures for compatibility
    specOut.data = SD.y;
    specOut.f = SD.f;
    specOut.peaksd = SP.y;
    specOut.peaksf = SP.f;
    specOut.peaksRes = SP.res;
    specOut.pSNR = SP.snr;
    specOut.dNsd = SP.dNsd;
    specOut.peaksFlag = SP.flag;
    specOutParams = Pfin;
    % convert some data to single to reduce storage
    specOut.data = single(specOut.data);
    specOut.peaksd = single(specOut.peaksd);
    specOut.peaksRes = single(specOut.peaksRes);
    specOut.pSNR = single(specOut.pSNR);
    % add extra info to parameters
    specOutParams.fileVersion = mfilename;
    % save
    save(fileName,'specOutParams','specOut');
    fprintf('done\n');
    clear specOut

    stitchfile_html{end+1} = opfname.stitchm;
    
    % save peak list as text file
    fprintf('Saving peak list text file...');
    if isempty(dir(fileList.peaksDir))
        fprintf('(making peaks dir)...');
        mkdir(fileList.peaksDir);
    end
    fileName = [fileList.peaksDir,'peakList_',fileList.setUID,'_',fileList.rawShort{fi},'.txt'];
    % save file
    fid = fopen(fileName,'wt');
    if fid==-1
        error('Cannot create file %s',fileName);
    end
    peaksmz = f2mz(SP.f, Pfin);
    if NOISE.NOISEFILT_ON
        % noise has been measured: output
        fprintf(fid,'M/Z\tINTENSITY\tSNR\tNON-NOISE_FLAG\n');
        for i=1:length(peaksmz)
           fprintf(fid,'%.7f\t%.6e\t%.2f\t%d\n',peaksmz(i),SP.y(i),SP.snr(i),SP.flag(i));
        end
    else
        % no noise measured
        fprintf(fid,'M/Z\tINTENSITY\n');
        for i=1:length(peaksmz)
            fprintf(fid,'%.7f\t%.6e\n',peaksmz(i),SP.y(i));
        end
    end
    fclose(fid);
    fprintf('done\n');

    peaklist_html{end+1} = ['peakList_',fileList.setUID,'_',fileList.rawShort{fi},'.txt'];
    
    
    % save internal calibration output
    if CAL.ON && fileList.calOn(fi)
        fprintf('Saving output from calibration...');
        fileName = [fileList.stitchDir,'calOutput_',fileList.setUID,'.txt'];
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
                fileList.rawShort{fi},...
                calibrants.ID{calOut(i).idx},...
                calibrants.mz(calOut(i).idx),...
                calOut(i).sw,...
                P(calOut(i).sw).mzStart,...
                P(calOut(i).sw).mzEnd,...
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
            fileList.rawShort{fi},...
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
    for i = 1:numSpectra
        fprintf(fid,'%s\t%i\t%2d\t%2d\t%2d\t%3d\t%3d\n',...
            fileList.rawShort{fi},...
            i,...
            P(i).dNnumpoints,...
            P(i).dNsd,...
            P(i).dNthresh,...
            P(i).dNrange(1),...
            P(i).dNrange(2));
    end
    fclose(fid);
    fprintf('done\n');
    
    

    
end

% ####################################################################
%% END MAIN LOOP
% ####################################################################

% ####################################################################
%% SAVE ALL 
% ####################################################################

%create html output for Galaxy
%write Stitch html for Galaxy
fid = fopen(html_outfile_stitch, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Process Transient Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>test</b></h2></div>');
fprintf(fid,'<hr></hr>')

for i=1:length(stitchfile_html)
  	html_code = ['<div style="padding: 3px"><b><a href="', stitchfile_html{i},'">',stitchfile_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'<p>');
fprintf(fid,'Note: ---');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);

%write peaklist html for Galaxy
fid = fopen(html_outfile_peaks, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Process Transient Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>test</b></h2></div>');
fprintf(fid,'<hr></hr>')
for i=1:length(peaklist_html)
  	html_code = ['<div style="padding: 3px"><b><a href="', peaklist_html{i},'">',peaklist_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'<p>');
fprintf(fid,'Note: ---');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);
fclose all;
exit;

%save parameters to text file
%fileName = message_file;
%fprintf('Saving parameters to file: %s\n',fileName);
%fid = fopen(fileName,'w');
%if ~fid, error('Cannot create message file'); end
%fprintf(fid,'FILE VERSION:\t%s',mfilename);
% debug
%fprintf(fid,'\r\n\r\nDEBUG:\r\n');
%fn = fieldnames(DEBUG);
%for i=1:length(fn)
%    val = DEBUG.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
% input files
%fprintf(fid,'\r\n\r\nFILES:\r\n');
%fn = fieldnames(FILES);
%for i=1:length(fn)
%    val = FILES.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
% noise filtering
%fprintf(fid,'\r\n\r\nNOISE FILTERING:\r\n');
%fn = fieldnames(NOISE);
%for i=1:length(fn)
%    val = NOISE.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
%fprintf(fid,'\r\n\r\nHIGHNOISEKNOWNREGIONS_MZ:\r\n');
%if ~NOISE.NOISEREMOVEALLKNOWN_ON
%    fprintf(fid,'\tSelected at runtime, defaults are:\r\n');
%else
%    for i=1:size(HIGHNOISEKNOWNREGIONS_MZ,1)
%        fprintf(fid,'\t%.5g\t%.5g\r\n',HIGHNOISEKNOWNREGIONS_MZ(i,:));
%    end
%end
% calibration
%fprintf(fid,'\r\n\r\nCAL:\r\n');
%fn = fieldnames(CAL);
%for i=1:length(fn)
%    val = CAL.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
% alignment
%fprintf(fid,'\r\n\r\nALIGN:\r\n');
%fn = fieldnames(ALIGN);
%for i=1:length(fn)
%    val = ALIGN.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else fc=fprintf(fid,'\t%s: %.10g \r\n',fn{i},val); end
%    if ~fc, error('Cannot write to message file'); end
%end
% stitch edge effects
%fprintf(fid,'\r\n\r\nSTITCH.NULL_REGION:\r\n');
%fn = fieldnames(STITCH.NULL_REGION);
%for i=1:length(fn)
%    val = STITCH.NULL_REGION.(fn{i});
%    if ischar(val), fc=fprintf(fid,'\t%s: %s \r\n',fn{i},val);
%    else
%        fc=fprintf(fid,'\r\n\t%s:',fn{i});
%        fc=fc & fprintf(fid,'\t%.10g',val);
%    end
%    if ~fc, error('Cannot write to message file'); end
%end
%fprintf('done\n');
%fclose(fid);
%fclose all;



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



function [spec,specParams] = GetRawProfileFS_v2(fileList,options,file,DISPLAY_HDRS,X_ZFILLS)
%returns data and frequency points (increasing values) from .raw profile file

%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHANGE HISTORY
%version 1.1: removed criteria for all A parameters to be the same: now just issues a warning
%version 1.2: updated for new fileList structure
% 25/Oct/07: Added message to press any key on warning
%version 1.3
% 13/Mar/08:    updated for new f2mz and mz2f functions, now delete Raw COM
%               server
%
% Current Version 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%global constants

if isempty(X_ZFILLS), error('X_ZFILLS not defined'); end

%intialise
spec = [];
specParams = [];

if ~isfield(options,'getSpec'), options.getSpec = 1; end
if ~isfield(options,'sumScans'), error('sumScans not specified'); end

%get Raw and Detector handles
%%%% hRaw = GetRawHandle(fileList.rawFull{file});
%%%% hDetector = GetDetectorHandle(hRaw);

if isunix
    system(['wine python ', which('MSFileReaderPy.py'), ' -i "' fileList.rawFull{file}, '" -o "', fullfile(fileList.datDir, 'temp.mat"')]);
else
    system([which('MSFileReaderPy.exe'), ' -i "' fileList.rawFull{file}, '" -o "', fullfile(fileList.datDir, 'temp.mat"')]);
end

matPy = load(fullfile(fileList.datDir, 'temp.mat'));

%get the list of scans for this spectrum
if isequal(fileList.scans{file},'all')
    %%%% scanList = hDetector.get('FirstSpectrum'):hDetector.get('LastSpectrum');
    scanList = matPy.temp.ScanList;
else
    scanList = fileList.scans{file};
end

%get handle to spectrum and header for each scan
%%%% hFilters = hDetector.get('Filters');
hFilters = matPy.temp.Filters;


%%%% filtersCount = hFilters.count;
filtersCount = matPy.temp.FiltersCount;

%if there is only one filter: easy
if filtersCount == 1
    filter = 1;
    %%%% hFilter = hFilters.Item(filter);
    hFilter = hFilters(1);
    %%%% if DISPLAY_HDRS, disp(['Filter: ', hFilter.Text]); end
    if DISPLAY_HDRS, disp(['Filter: ', hFilter]); end
    filterOfScanIndex = repmat(filter,1,length(scanList));
%otherwise need to find the filter index for each scan
else
    for i=1:length(scanList)
        %%%% hFilter = hFilters.ScanNumber(scanList(i));
        hFilter = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(scanList(i)) ,'.Filter'));
        targetFilter = hFilter;
        %%%% targetFilter = hFilter.Text;
        %%%% if DISPLAY_HDRS, disp(['Filter: ', hFilter.Text,', scan: ',num2str(i)]); end
        if DISPLAY_HDRS, disp(['Filter: ', hFilter,', scan: ',num2str(i)]); end
        
        found = 0;
        for filterIndex = 1:filtersCount
            %%%% hFilter = hFilters.Item(filterIndex);
            hFilter = hFilters(filterIndex);
            %%%% if isequal(hFilter.Text,targetFilter), found=1; break; end
            if isequal(char(hFilter),targetFilter), found=1; break; end
        end
        if ~found, error('Cannot find filter in GetSpectraHandle.m'); end
        filterOfScanIndex(i) = filterIndex;
% THE FOLLOWING TAKES UP TOO MUCH MEMORY HERE
%         hSpectra = hDetector.get('Spectra', filterIndex);
%         hSpectrum(i) = hSpectra.ScanNumber(scanList(i));
%         hSpectrumHeader(i) = hSpectrum(i).get('Header');  %get header for parameters
    end
end

%the list and count of filters actually used is now updated
filtersUsed = unique(filterOfScanIndex);
filtersCount = length(filtersUsed);

%get the centroid (labeldata)
% if options.sumScans
%     %concat the centroids from each scan
%     spectrum.centroidsMZ = [];
%     spectrum.centroidsData = [];
%     for i=1:numberOfScans
%         temp = hSpectrum(i).get('LabelData');
%         spectrum.centroidsMZ = [spectrum.centroidsMZ temp(1,:)];
%         spectrum.centroidsData = [spectrum.centroidsData temp(2,:)];
%     end
%     [spectrum.centroidsMZ,idx] = sort(spectrum.centroidsMZ);
%     spectrum.centroidsData = spectrum.centroidsData(idx);
% elseif options.stitch
    spectrum.centroidsMZ = [];
    spectrum.centroidsData = [];
% else
%     for i=1:numberOfScans
%         temp = hSpectrum(i).get('LabelData');
%         spectrum(i).centroidsMZ = temp(1,:);
%         spectrum(i).centroidsData = temp(2,:);
%     end
% end

%get the spectra
if DISPLAY_HDRS, disp(['File: ',fileList.rawShort{file},'...']); end
tic; lastT=0;

%loop through each filter (these will form the spectra)
spectraScans = {};
spectraFilter = {};

%options.sumScans = 1;
if options.sumScans
    
    for k = 1:filtersCount

        filter = filtersUsed(k);
        scansThisFilter = [];
        d = [];

        %first get a list of scans from this filter
        idx = find(filterOfScanIndex==filtersUsed(k));
        filterScans = scanList(idx);
        scanIT = [];
        scanTIC = [];

        %get the handle to the spectra for this filter
        %%%% hSpectra = hDetector.get('Spectra', filter);

        %loop through each scan in this filter
        for i=1:length(filterScans)

            %get the spectrum for this scan
            %%%% hSpectrum = hSpectra.ScanNumber(filterScans(i));
            %%%% hSpectrumHeader = hSpectrum.get('Header');

            %first check for varying conversion parameters by extracting the trailer extra info
            %%%% hTrailerExtras = hDetector.TrailerExtras();
            %%%% hTrailerExtra = hTrailerExtras.Item(filterScans(i));
            %%%% trailerExtraData = hTrailerExtra.get('Data');
            %%%% newParams.A = str2num(trailerExtraData{2,25})*1e3;
            newParams.A = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.A'))*1e3;
            if i==1, fParams.A = newParams.A;
            elseif newParams.A~=fParams.A
                warning(['Parameter A is changing in ',fileList.rawShort{file},' - included, but check scan ',num2str(filterScans(i))]);
            end
            
            %%%% newParams.B = str2num(trailerExtraData{2,26})*1e6;
            newParams.B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
            if i==1, fParams.B = newParams.B;
            elseif newParams.B~=fParams.B
                warning(['Parameter B is changing in ',fileList.rawShort{file},' - included, but check scan ',num2str(filterScans(i))]);
            end

            
            %%%% scanTIC = [scanTIC hSpectrumHeader.TIC];
            scanTIC = [scanTIC eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.TIC'))];

            %%%% scanIT = [scanIT str2num(trailerExtraData{2,3})];
            scanIT = [scanIT eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'))];
            % eval(strcat('matPy.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'));
            
            %%%% if str2num(trailerExtraData{2,3}) >= 750
            if eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT')) >= 750
                warning(['Underfill in ',fileList.rawShort{file},' - NOT!!! ignoring scan ',num2str(filterScans(i))]);
                %continue;
            end
            
            if options.getSpec
                %get spectrum
                %%%% tempSpectrum = hSpectrum.get('Data');

                %can store the data values as singles to save space - can't do the same with m/z values, need the precision

                %%%% scanmz = tempSpectrum(1,:);
                scanmz = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scanmz'));
                %%%% scand = tempSpectrum(2,:);
                scand = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scand'));
                
                % to extract the noise data, one would used
                % hSpectrum.get('NoiseData');
                % which returns 4 rows: m/z, noise and baseline.
                
                %if the current A or B parameters have changed (are not the selected global parameters), normalise the m/z values to the global A and B parameter
                if newParams.B~=fParams.B || newParams.A~=fParams.A
                    scanmz = mz2f(scanmz,newParams);
                    scanmz = f2mz(scanmz,fParams);
                end
                
                %store scan data
                if isempty(d)
                    d = scand;
                    mz = scanmz;
                else    %append data
                    newd = [];
                    for j=1:size(d,1)
                        %in the case that there are multiple scans, we need to keep the
                        %number of data point entries the same, such that all d's can
                        %use the same 'mz'
                        [newd(j,:),scand_temp,newmz] = FillOut(d(j,:),mz,scand,scanmz);
                    end
                    d = sum([newd;scand_temp],1);   %sum here to increase efficiency
                    mz = newmz;
                end
            end

            scansThisFilter = [scansThisFilter filterScans(i)];

        end

        spectraScans{k} = scansThisFilter;
        spectraFilter{k} = filter;

        %average the parameters if summing spectra over each filter: save each value otherwise
        IT(k) = mean(scanIT);
        TIC(k) = mean(scanTIC);
        
        %%%% mzStart(k) = hSpectrumHeader.LowMass;
        mzStart(k) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.LowMass'));
		
        %%%% mzEnd(k) = hSpectrumHeader.HighMass;
        mzEnd(k) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.HighMass'));
        
        params(k) = fParams;
%         else
%             IT((end+1):specCount) = scanIT;
%             TIC((end+1):specCount) = scanTIC;
%             mzStart((end+1):specCount) = hSpectrumHeader.LowMass;
%             mzEnd((end+1):specCount) = hSpectrumHeader.HighMass;
%             params((end+1):specCount).A = fParams.A;
%             params((end+1):specCount).B = fParams.B;
%         end

        if options.getSpec
            %average spectra abundance
            data = d ./ length(scansThisFilter);
            %convert to frequency (increasing values); decalibrate
            spectrum(k).f = mz2f(mz,params(k));
            spectrum(k).f = fliplr(spectrum(k).f);
            spectrum(k).data = fliplr(data);
        end

        %progress counter
        progress = k/filtersCount*100;
        t=toc;
        if t<5, lastT = floor(t/10);
        elseif t>5 && floor(t/10)>lastT
            disp([num2str(progress),'%']);lastT=floor(t/10);drawnow;
        end
        
    end
    numberOfSpectra = filtersCount;

else    %not summing scans

    scanCount = 0;;
    for k = 1:filtersCount

        filter = filtersUsed(k);

        %first get a list of scans from this filter
        idx = find(filterOfScanIndex==filtersUsed(k));
        filterScans = scanList(idx);

        %get the handle to the spectra for this filter
        %%%% hSpectra = hDetector.get('Spectra', filter);

        %loop through each scan in this filter
        for i=1:length(filterScans)

            %get the spectrum for this scan
            %%%% hSpectrum = hSpectra.ScanNumber(filterScans(i));
            %%%% hSpectrumHeader = hSpectrum.get('Header');
            
            scanCount = scanCount + 1;
            spectraScans{scanCount} = filterScans(i);
            spectraFilter{scanCount} = filter;

            %%%% hTrailerExtras = hDetector.TrailerExtras();
            %%%% hTrailerExtra = hTrailerExtras.Item(filterScans(i));
            %%%% trailerExtraData = hTrailerExtra.get('Data');
            
            %%%% TIC(scanCount) = hSpectrumHeader.TIC;
            TIC(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.TIC'));
            %%%% IT(scanCount) = str2num(trailerExtraData{2,3});
            IT(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.IT'));
            %%%% mzStart(scanCount) = hSpectrumHeader.LowMass;
            mzStart(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.LowMass'));
            %%%% mzEnd(scanCount) = hSpectrumHeader.HighMass;
            mzEnd(scanCount) = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.HighMass'));
            %%%% params(scanCount).A = str2num(trailerExtraData{2,25})*1e3;
            params(scanCount).A = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.A'))*1e3;
            %%%% params(scanCount).B = str2num(trailerExtraData{2,26})*1e6;
            params(scanCount).B = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.B'))*1e6;
            
            if options.getSpec
                %get spectrum
                %%%% tempSpectrum = hSpectrum.get('Data');

                %can store the data values as singles to save space - can't do the same with m/z values, need the precision
                %%%% scanmz = tempSpectrum(1,:);
                scanmz = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scanmz'));
                %%%% scand = tempSpectrum(2,:);
                scand = eval(strcat('matPy.temp.ScanInfo.Scan', num2str(filterScans(i)) ,'.scand'));

                %progress counter
                progress = i/length(filterScans)*100;
                t=toc;
                if t<5, lastT = floor(t/10);
                elseif t>5 && floor(t/10)>lastT
                    disp([num2str(progress),'%']);lastT=floor(t/10);drawnow;
                end

                %convert to frequency (increasing values); decalibrate and store
                spectrum(i).f = mz2f(scanmz,params(i));
                spectrum(i).f = fliplr(spectrum(i).f);
                spectrum(i).data = fliplr(scand);
            end
        end
    end
    numberOfSpectra = scanCount;
end
    
%now pad out (uncompress) the data and frequency values
%theoretical delta f is BW/n, but we don't know BW/n
if options.getSpec
    for i=1:numberOfSpectra

        %first estimate the best frequency step to use (there are slight variations
        %in the actual frequency step probably due to representation and rounding errors.
        deltaf = median(spectrum(i).f(2:end)-spectrum(i).f(1:end-1));
        T(i) = 1/deltaf/2^X_ZFILLS;
        if round(T(i)*1000) ~= 768, warning(['T not 0.768 in raw file, is ',num2str(T(i))]); keyboard; end;

        %each data point will have a new index relative to the first frequency
        newIdx = (spectrum(i).f-spectrum(i).f(1))./deltaf;  %index relative to first freq. val
        if find(abs(round(newIdx)-newIdx)>0.3), warning(['Possible error finding data indices in GetRawProfileFS.m! Press any key to continue...']); pause; end
        newIdx = round(newIdx) + 1;

        %spread out the data into the new data matrix, at the same time padding the gaps with zero
        spectrum(i).newData(newIdx) = spectrum(i).data;

        %the new frequency values are the indices multiplied by the frequency step, plus the starting frequency
        spectrum(i).newf = ((0:length(spectrum(i).newData)-1).*deltaf) + spectrum(i).f(1);
    end
end

%return values of interest
numberOfSpectra;
for i=1:numberOfSpectra
    if options.getSpec
        spec(i).y = spectrum(i).newData;
        spec(i).f = spectrum(i).newf;
        specParams(i).T = T(i);
    end
%     spec(i).xCentroidsMZ = spectrum(i).centroidsMZ;
%     spec(i).xCentroidsData = spectrum(i).centroidsData;
    specParams(i).mzStart = mzStart(i);
    specParams(i).mzEnd = mzEnd(i);
    specParams(i).TIC = TIC(i);
    specParams(i).A = params(i).A;
    specParams(i).B = params(i).B;
    specParams(i).C = 0;
    specParams(i).zFills = X_ZFILLS;
    specParams(i).IT = IT(i);
    specParams(i).filterIndex = spectraFilter{i};
    specParams(i).scans = spectraScans{i};
end

% release all interfaces derived from the Raw COM server
% and then delete the server itself.
%%%% hRaw.delete;


function [hRaw] = GetRawHandle(fileName)
%returns a new RAW handle given a filename

% MODIFICATION HISTORY
% 14/Mar/08         Added error code output.

try
    hRaw = actxserver('XRaw.XRaw');
    hRaw.Open(fileName);
catch
    errCode = lasterr;
    errCode = hex2dec(errCode(end-3:end));
    fprintf('Error loading file: %s\n',fileName);
    switch errCode
        case 10001
            error('File not found.');
        case 10002
            error('Internal error: wrong dll version? Try running reg.m');
        case 10006
            error('Cannot access information in the specified file type.');
        case 10007
            error('The file is already open.');
        case 10008
            error('No FileName supplied, or the name of the file was not set up in another method.');
        case 10009
            error('An invalid parameter value has been supplied.');
        case 10012
            error('The file is of an incorrect type.');
        case 10014
            error('There is no valid data available.');
        case 10017
            error('The object is not valid.');
        otherwise
            rethrow(lasterror);
    end
end


function [hDetector] = GetDetectorHandle(hRaw)
%gets the detector handle for the raw handle

%constants
MAX_DETECTORTYPES = 5;

%list all possible Detector Types
hRawInfo = hRaw.get('RawInfo');
DetectorCount = hRawInfo.get('DetectorCount');
if DetectorCount > 1
    disp(['Number of detectors in .raw: ', num2str(DetectorCount)]);
    %DetectorsOfTypeCount
    %there are 6 possibilities to detectortype. seems
    %XNo-Device is not an option
    for i = 0:(MAX_DETECTORTYPES-1)
        DetectorsOfTypeCount(i+1) = hRawInfo.get('DetectorsOfTypeCount',i);
    end
    disp('Detectors Of Type Count: '); disp(DetectorsOfTypeCount);
    detectorTypes = hRaw.set('Detector');
    disp('Detector Types: '); disp(detectorTypes);
    reply = input('Which Detector Type (0 to 5)? [0]:');
    if isempty(reply)
        reply = 0;
    end
else
    reply = 0;
end

% XDetectorRead based on detector type choice
try
    hDetector = hRaw.get('Detector',reply,1);
catch
    Xerrors(lasterr);
    error('Error in GetDetectorHandle');
end


function f = mz2f(mz, P, y)
%converts from mz, returns frequency
%P.A, P.B and P.C contain calibration parameters
%if C is supplied, it is expected to apply to the form of calibration equation derived by Masselon et al. 2002, ref 63:
% m/z = A./f + B./f.^2 + C*I./f.^2, where I is the intensity of each peak

numTerms = 2;   % number of parameters in calibration equation

if nargin == 1
    % no P parameter supplied
    disp('WARNING: Insufficient parameters specified, using nominal values');
    numTerms = 2;
    A = 107368.5e3;
    B = -750.461e6;
elseif isempty(P)
    % P is empty: use defaults
    numTerms = 2;
    A = 107368.5e3;
    B = -750.461e6;
elseif ~isfield(P,'C') || P.C==0
    % no (or zero) C parameter: use two-term calibration
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
        f = (A + sqrt( A^2+mz.*4*B ))./(2.*mz);
    case 3
        %solve using symbolic
        syms As Bs Cs Is Fs mzs
        mz = -mzs + As./fs + Bs./fs.^2 + Cs*ys./fs.^2;
        f = solve(mz,fs);
        mzL = length(mz);
        if mzL > 1
            A = A.*ones(1,mzL);
            B = B.*ones(1,mzL);
            C = C.*ones(1,mzL);
        end
        f = (subs(f,{As,Bs,Cs,ys,mzs},{A,B,C,y,mz}));
        f = max(real(f));
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



function [d1,d2,mz2] = FillOut(d1,mz1,d2,mz2)
%make both inputs with same frequency points, extra data points are set to
%zero

%merge and sort
mz1 = [mz1 mz2];
[mz1 idx]=sortrows(mz1.');
mz1 = mz1.';

%sort data same way
d2 = [zeros(1,length(d1)) d2];
d1(length(mz1)) = 0;

d1 = d1(idx);
d2 = d2(idx);

%remove any duplicates
idx = find((mz1(1:(end-1))-mz1(2:end))==0);
if ~isempty(idx)
    mz1(idx) = [];
    d1(idx+1) = [];
    d2(idx) = [];
end

mz2 = mz1;



function [rms,fitSigma,fRange, numPoints] = GetNoiseLevel_0_4(data,f,exregions)
% This function estimates the noise level of data, assuming the noise
% follows a Rayleigh distribution and that the data is predominantly noise.
% A section of test data is iteratively picked and fitted to the Rayleigh
% distribution.  The threshold (on the PDF) corresponding to the presence
% of a single data value above that threshold is determined, and this is
% compared to the actual maximum value in the test range selected.  If the
% actual maximum is higher, the region contains signal and the next region
% is selected.
% The test region is first moved up the data by a jump of delta, if no
% noise region is found and the end of the data is reached, the test region
% is reduced in size by rangeScale, and the process iterates.

% INPUTS
% data, f: mass spectrum data points
% exregions: frequency ranges to exclude from search ([range1;range2;...])

% OUTPUTS
% rms: rms noise level (of actual data)
% fitSigma: Rayleigh parameter (of fitted data)
% fRange: frequency range [from,to] over which noise is measured

% version 0.1, 20/Sep/07 - from GetNoiseThreshold_0_1, this function now simply returns the noise level
%                          (ie does not calculated an adjusted threshold - this moved to Stitch).
%              16/Jan/08 - displays m/z range over which noise calculated.
% version 0.2  06/Feb/08 - returns freq range over which noise calculated.
%              25/Feb/08 - minor changes, added parameters
% version 0.3  28/Feb/08 - catches case where noise region overlaps with
%                           region of "known" noise and searches for
%                           alternative
%              10/Mar/08 - changed delta to "delta = round(numPoints/10);"
%                           from "delta = round(numPoints/2);"

MAXPOINTS = 10000;   % maximum (starting) number of data points for parameter estimation
MINPOINTS = 50;    % minimum number of data points for parameter estimation

numPoints = min([MAXPOINTS length(data)]); % starting number of data points for noise estimation
delta = round(numPoints/100); % jump in index
rangeScale = 1.01;   % divisor for reducing the length of test range

%find a suitable range data with no real peaks
done = 0;   % flag
iStart = 1; % starting index
while ~done
    %select data points
    testData = data(iStart:iStart+numPoints-1);
    testf = f(iStart:iStart+numPoints-1);

    %fit
    ray = fit_ML_rayleigh(testData); %pause
    %disp(ray.RMS)
    maxVal = max(testData);

    %max limit
    %P(X<=x) = 1-exp(-x^2 / 2sigma^2)
    %solving gives: x = sigma . sqrt(-2log(1-P))
    %1-P = 1/numPoints, ie on average one data point will exceed the limit
    sigma = sqrt(ray.s);
    maxLim = sigma*sqrt(-2 * log(1/numPoints));

    %move test region along by delta
    iStart = iStart + delta;
    
    %check if noise region includes exregions
    if find(min(exregions,[],2) < max(testf) & max(exregions,[],2) > min(testf))
        inEx = 1;
    else
        inEx = 0;
    end

    if (maxVal <= maxLim) & ~inEx
        % found region of noise
        fprintf('Found signal-free region over %d data points (%.2f-%.2fm/z), sigma = %0.2f\n',numPoints,f2mz(testf(end),[],[]),f2mz(testf(1),[],[]),sigma);
        % debug output
        fprintf('Sigma from mean = %0.2f, rms error = %.2f\n',mean(testData)/sqrt(pi/2),ray.RMS);
%         figure;h = gca;hist(h,testData,100);
%         plot(f2mz(f,[],[]),data);axis([f2mz([max(testf) min(testf)],[],[]) 0 max(testData)]);
%         close;
        % end debug output
        done = 1;
    elseif (iStart + numPoints > length(data))
        % got to end of data: decrease test region length and reset to
        % start of data
        iStart = 1;
        numPoints = round(numPoints/rangeScale);
        delta = round(numPoints/2);
        if numPoints < MINPOINTS
            plot(f2mz(f,[],[]),data);
            error('Can''t find long enough region with no signal: check spectrum');
        end
    end    
end

rms = norm(testData)/sqrt(numPoints);
fitSigma = sigma;
fRange = [min(testf) max(testf)];


function result = fit_ML_rayleigh( x,hAx )
% fit_ML_rayleigh - Maximum Likelihood fit of the rayleigh distribution of i.i.d. samples!.
%                  Given the samples of a rayleigh distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(r)=r*exp(-r^2/(2*s))/s
%    with parameter: s
%
% format:   result = fit_ML_rayleigh( x,hAx )
%
% input:    x   - vector, samples with rayleigh distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      s   - fitted parameter
%                      CRB - Cram?r-Rao Bound for the estimator value
%                      RMS - RMS error of the estimation 
%                      type- 'ML'
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The rayleigh destribution is given by:
%
%    p(x,s) = x / s * exp(-x^2/(2*s))
%
%    where x are the samples which distribute by the function p(x,s)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    f(Xn,s)   = Xn / s * exp(-Xn^2/(2*s))
%    L(s)      = f(X,s) = product_by_n( f(Xn,s) )
%              = PI(Xn) * (s^(-N)) * exp( -sum(Xn^2)/(2*s) )
%    log(L(s)) = sum(log(Xn)) - N*log(s) - sum(Xn^2)/(2*s)
%     
%    The maximum likelihood point is found by the derivative of log(L(s)) with respect to "s":
%
%    diff(log(L(s))) = -N/s + sum(Xn^2)/(2*s^2) = (N/s^2) * ( sum(Xn^2)/(2*N) - s ) = 
%                    = J(s) * (s_estimation - s)  
%
%    Therefore, the (efficient) estimator is given by:
%
%               s = sum( Xn^2 ) / (2 * N)
%
%    The Cram?r-Rao Bound for this estimation is:
%
%               VAR( s ) = 1/J(s) = (s^2)/N
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

if (nargin<1)
    error( 'fit_ML_rayleigh - insufficient input arguments' );
end

% Estimation
% =============
x       = real(x(:));                 % should be column vectors !
N       = length(x);
s       = sum(x.^2)/(2*N);
CRB     = (s^2)/N;
[n,x_c] = hist( x,100 );
n       = n / sum(n*abs(x_c(2)-x_c(1)));
y       = x_c.*exp(-x_c.^2/(2*s))/s;
RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );

% finish summarizing results
% ============================
result = struct( 's',s,'CRB',CRB,'RMS',RMS,'type','ML' );

% plot distribution if asked for
% ===============================
if (nargin>1)
    xspan = linspace(min(x),max(x),100);
    if ishandle( hAx )
        plot_rayleigh( xspan,result,hAx,3 );
    else
        figure;
        plot_rayleigh( xspan,result,gca,3 );
    end
end



function [vals]=maxima(input_vec, number);

%
% maxima - return the locations of a given number of maxima 
%
% [output] = maxima(input, number of maxima)
%
% inputs : input_vec = a vector containing the input data
%          number = (obviously) the number of maxima to find
% 
% output : vals = vector containing the locations of the maxima found
%
% GPP, Uni of B'ham, 2001

%Do a quick check on the input incase the user's a fool
if ((size(input_vec, 1) ~= 1) & (size(input_vec, 2) ~= 1) | (ndims(input_vec) > 2))
    disp('Stop being a muppet. The input cannot be a matrix');
    return;
end

%Turn the input into a row-vector
if (size(input_vec, 1) > 1)
    input_vec = input_vec.';
end

%Firstly, find all the (local) maxima in the input
deriv = input_vec(2: size(input_vec,2)) - input_vec(1:size(input_vec,2)-1);
deriv_sq = deriv(2: size(deriv,2)) - deriv(1:size(deriv,2)-1);

mult = deriv(2: size(deriv,2)) .* deriv(1:size(deriv,2)-1);
turning_points = find(mult <= 0);

maxima = turning_points(find(deriv_sq(turning_points) < 0));

%Now, check if we have enough maxima - if not, just output those that we do have:
if (size(maxima,2) <= number)
    vals = maxima;
    return;
end

%Next, apply a value row above the maxima & sort the columns (ascending order)
maxima = sortrows([input_vec(maxima+1);maxima].').';

%Finally, output the locations of the biggest 'number' of rows
vals = fliplr(maxima(2,size(maxima,2)-number + 1:size(maxima,2)));

function [peaksy, peaksh, peaksf, peaksr] = PeakDetectArea_0_3(y,f,mode,thresh)
%find the peaks in data
% 
%output is
%-peaks .y (peak area)
%       .h (height)
%       .f (centre frequency)
%       .r (peak resolution)
% 
%inputs are:
%-f         frequency axis in Hz, increasing values
%-y         corresponding data points
%-mode      peak detect mode
%-thresh    minimum peak height - below this, peaks will not be returned

%CHANGES
% 22/Jun/07 Filter out peaks with <3 non-zero data points, area
%           determination not working for these (insufficient data for
%           determining peak area and center anyway)
% 01/Aug/07 (KCEinterp) include data points either side of maxima even if below noise
% 6/Sep/07  v0.1    Now also calculates the FWHM for peaks above the noise
%           floor and code optimised
% 6/Feb/08  v0.2    Added possibility to not apply peak threshold
%                   Returns maximum (existing) data point for peaks
% 14/Feb/08         thresh is now the maximum allowed data point (not peak
%                   height)
% 01/Apr/08         Updated centroid option. Removed QUIET. Removed
%                   redundant peak detection options
% 08/Apr/08 v0.3    Returns peak height instead of max. data point, and slightly optimised

switch mode
    case 'centroid'
        %centroid - computers in MS p35
        % quick: doesn't really work but gives an approximate answer
        [peaksy, peaksh, peaksf, peaksr] = Centroid(y,f,thresh);
    case 'KCEinterp'   %ref 56
        KCexp = 5.5;
        % other values of the KC exponent can be used for different line shapes - see ref 56
        [peaksy, peaksh, peaksf, peaksr] = KCeInterpolate(y,f,KCexp,thresh);
    otherwise
        error('Mode choice not recognised!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peaksy, peaksh, peaksf, peaksr] = KCeInterpolate(y,f,e,thresh)
% 
% INPUT
% y, f:     spectral data points
% e:        KC exponent
% thresh:   min peak height to return
% 
% OUTPUT
% see above

% find indices of maxima
%  limitations:
%  maxima at postion (1) or (end) are not included
%  exact adjacent data counted as 2 peaks eg [0 1 1 1 0] -> peaks at [2 4]
peaksi = maxima(y,length(y)) + 1;
N = length(peaksi);

peaksy = zeros(1,N);
peaksh = zeros(1,N);
peaksf = zeros(1,N);
peaksr = zeros(1,N);

% loop through each peak
fprintf('Progress:  0%%');
tic;
n = 0;
for i=1:N
    %fit: 3 data points
    j = (peaksi(i)-1):(peaksi(i)+1);
    X = [];
    dfoff = f(j(2));   %avoid rounding errors in inversion by using relative values
    foff = f(j) - dfoff;
    yy = y(j);
    X(:,1) = (foff.^2).';
    X(:,2) = foff.';
    X(:,3) = ones(3,1);
    Y = (yy.^(1/e)).';
%         Z = X\Y    %Gaussian elimination
    [U,S,V] = svd(X);
    Z = (V*inv(S)*U')*Y;
    a = Z(1);
    b = Z(2);
    c = Z(3);
    %fitted peak centre
    x0 = foff(2) - abs(foff(2)-foff(1))/2 * (yy(3)-yy(1)) / (yy(1) - 2*yy(2) + yy(3));
    %find peak height
    y0 = (a*x0.^2 + b*x0 + c)^e;
    
    % is peak wanted?
    if y0 > thresh
        % yes
        n = n+1;
        peaksf(n) = x0 + dfoff;
        peaksh(n) = y0;
        % calculate the resolution FWHM
        q = sqrt(b^2-4*a*(c-exp(log(y0/2)/e)));
        peaksr(n) = f2mz(x0+dfoff,[]) / abs(f2mz((-b+q)/(2*a)+dfoff,[]) - f2mz((-b-q)/(2*a)+dfoff,[]));
        % accurately integrate using KCe interpolation
        q = 4*a*c - b^2;
        qq = sqrt(-q);
        zeroX(1) = (-b+qq)/(2*a);
        zeroX(2) = (-b-qq)/(2*a);
        if length(zeroX)~=2, disp(zeroX); error('Check zero-crossings'); end
        %sort limits
        if zeroX(1)>zeroX(2), zeroX = fliplr(zeroX); end
        %integrate using maple - slow
%         peaksd(n) = (double(int(ys,xs,zeroX(1)-foff,zeroX(2)-foff)));
        %integrate using indefinite solution - as per paper- fast
        z = e - 0.5;
        k = 4*a/q;          %NOTE MISPRINT IN PAPER!!!! ARRRGGHGHHH!!!  Took me hours to find that!
        ya = XRootXIntegral(z,k,a,b,c,zeroX(2),q) - XRootXIntegral(z,k,a,b,c,zeroX(1),q);
        %may be small imaginary component due to rounding errors (allow 0.1%)
        if abs(imag(ya)/real(ya)) > 0.001, error('Check imag components of peak area'); end
        peaksy(n) = real(ya);
    end

    % progress
    if floor(toc) || (i==N) % update every second
        fprintf('\b\b\b');
        fprintf('%2.0f%%',floor(i/N*100));
        tic;
    end
end
fprintf('\n');
% resize matrices
peaksy = peaksy(1:n);
peaksh = peaksh(1:n);
peaksf = peaksf(1:n);
peaksr = peaksr(1:n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plotting codefrom KCeInterpolate
        %plot curve fitting to
%             j = (peaksn(i)-10):(peaksn(i)+10);
%             j = j(find(j>0));
%             figure;plot(f2mz(f(j)),data(j),'k'); hold all;
        %plot fitted peak cetnre - KCe
%             stem (f2mz(x0),y0,':r');
        %plot fitted peak cetnre - GIalpha
%         j = (peaksn(i)-1):(peaksn(i)+1);
%         r = (data(j(1))/data(j(2)))^2;
%         s = (data(j(3))/data(j(2)))^2;
%         k = f(j(2));
%         alpha = 2.0;
%         x0GI = k - ((alpha+1)*(r-s)/(2*(alpha+r+s-r*s*(alpha+2))));
%         stem (x0GI,data(j(2)),':g');
        %plot fitted function
%             j = (peaksn(i)-20):(peaksn(i)+20);
%             j = j(find(j>0));                 %incase we've gone beyond the beginning of the spectrum
%             X = linspace(max(f(j)),min(f(j)));
%             Y = (a*(X-foff).^2 + b*(X-foff) + c);
%             Y(find(Y<0)) = 0;
%             Y = Y.^e;
%             plot(f2mz(X),Y,':r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [myInt] = XRootXIntegral(z,k,a,b,c,x,q)
% determines the integral of the function R^z.sqrt(R)
% at the limit x 
% where R = a*x^2 + b*x + c
% and q = 4*a*c - b^2
% and k = 4*a/q
% 
% NOTE: code optimised for x specified at y=0 values: will not work for
% other values of x!!

if a>0
    inta = 1/sqrt(a)*asinh((2*a*x+b)/sqrt(q));
else
    inta = -1/sqrt(-a)*asin((2*a*x+b)/sqrt(-q));
end
R = (a*x^2 + b*x + c);
%     suma = sum(factorial([0:z]) .* factorial([1:z+1]) .* (4*k*R).^([0:z]) ./ factorial(2*[0:z]+2));
% note that since we are looking for the zero-crossing, R will always
% be tending ~0, and so suma ~ 0.5
suma = 0.5;
myInt = factorial(2*z + 2) / ((factorial(z+1))^2 * (4*k)^(z+1)) * ( k*(2*a*x + b)*sqrt(R) / a * suma + inta );
%     myInt = factorial(2*z + 2) / ((factorial(z+1))^2 * (4*k)^(z+1)) * inta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peaksy, peaksh, peaksf, peaksr] = Centroid(y,f,thresh)
%use the centroiding to return peaks
%uses only the maximum point and the two data points either side.

% indices of maxima
peaksi = maxima(y,length(y)) + 1;

% maximum data point for each peak
peaksm = y(peaksi);

% frequencies
peaksf=y(peaksi-1).*f(peaksi-1) + y(peaksi).*f(peaksi) + y(peaksi+1).*f(peaksi+1);
peaksf=peaksf./(y(peaksi-1) + y(peaksi) + y(peaksi+1));
%peaksf=data(peaksn-2).*f(peaksn-2) + data(peaksn-1).*f(peaksn-1) + data(peaksn).*f(peaksn) + data(peaksn+1).*f(peaksn+1) + data(peaksn+2).*f(peaksn+2);
%peaksf=peaksf./(data(peaksn-2) + data(peaksn-1) + data(peaksn) + data(peaksn+1) + data(peaksn+2));

% peak height (actually peak data maxima instead)
peaksh = peaksm;

% peak area (actually peak data maxima instead)
peaksy = peaksm;

% filter if required
if thresh
    idx = peaksh>thresh;
    peaksy = peaksy(idx);
    peaksf = peaksf(idx);
    peaksh = peaksh(idx);
end

% peak resolution: set to zero
peaksr = zeros(size(peaksy));

fprintf('\n');


function [idxTarget,idxIn] = FindClosest_0_3(target,inputList,maxErr,errType)
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


function [idxTarget,idxIn] = FindLargest_1_0(target,inputList,maxErr,errType, intensities)
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
%070912 (becomes FindLargest_1_0) - intensities is new input. Needs errtype
%                   to be ppm and maxErr to be defined. Checks all peaks in range and selects
%                   largest

%check for empty input
if (isempty(inputList) | isempty(target)); warning('Empty input'); idxTarget=[]; idxIn=[]; return; end;


%{
%sort inputList
[inputList,idxSortInput] = sort(inputList.');
inputList = inputList.';

%sort target
[target,idxSortTarget] = sort(target.');
target = target.';
%}

for i = 1:length(inputList)
    MZth = inputList(i);
    MZrange = (MZth/1e6)*maxErr;
    MZmin = MZth-MZrange;
    MZmax = MZth+MZrange;
    index_peaks = intersect(find(target>=MZmin),find(target<=MZmax));
    intensities_tmp= intensities(index_peaks);
    [int,ind] = max(intensities_tmp);
    if ~isempty(ind)
        cal_ind(i) = index_peaks(ind);
        cal_mz(i) = target(cal_ind(i));
    else
        cal_ind(i)=1;
    end
end

idxTarget= cal_ind;
idxIn = 1:length(inputList);

%{



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
%}



function [P] = Calibrate_1_0(SP, P, calmz, mode, weighted)
%Determines new calibration parameters given observed peak frequencies (SP) and a list of
%known masses ie cal points (calmz) in m/z, and initial A and B parameters

% INPUTS
% SP is a structure with fields .y and .f corresponding to the peaks in the
%   spectrum identified as calibrants, SP sorted by increasing f.
% P is the parameters for SP
% calmz is the calibrant exact masses, must correspond to order of peaks in
%   SP
% mode is the calibration mode (b or ab or abc)
% weighted is 1: weight calibration to more intense peaks, 0: no weighting

% OUTPUTS
% P is updated with the new (recalibrated) calibration parameters

% UPDATES
% 18/6/07 (v0.1): FindClosest now v0.1
% 19/6/07       : Max peak distance now parameter
% 12/10/07      : Now also returns "data", the abundance of the calibrant
%                   in the spectrum
% 26/02/08 (v1.0):  Re-structured, removed output, input peaks are now
%                   already identified with calibrants

%check calmz in single column
if (size(calmz,2)>size(calmz,1)), calmz=calmz.'; end

%weighting
if weighted, calmzw = calmz.*(SP.y).'; end

%construct F for mz = A/f + B/f^2
f = SP.f.';
y = SP.y.';
A = P.A;
B = P.B;
C = P.C;
switch mode
    case 'b'
        F = 1./f.^2;
        if weighted
            F = F.*y;
            Xh = F\(calmzw - A*1./f.*y);
        else
            Xh = F\(calmz - A*1./f);
        end
        newA = A;
        newB = Xh;
        if C, error('Non-zero C value in b-mode recalibration?'); end
        newC = 0;
    case 'ab'
        F = [1./f 1./f.^2];
        if weighted
            F = F.*[y y];
            Xh = F\calmzw;
        else
            Xh = F\calmz;
        end
        newA = Xh(1);
        newB = Xh(2);
        if C, error('Non-zero C value in ab-mode recalibration?'); end        
        newC = 0;
    case 'abc'  %see ref 63, Masselon et al., 2002
        F = [1./f 1./f.^2 y./f.^2];
        if weighted
            F = F.*[y y y];
            Xh = F\calmzw;
        else
            Xh = F\calmz;
        end
        newA = Xh(1);
        newB = Xh(2);
        newC = Xh(3);
end

P.A = newA;
P.B = newB;
P.C = newC;



function resid = AlignMinFunc_0_1(rL,sL,rR,sR,params)
% this function is the minimisation function for the alignment of
% overlapping SIM windows. The function is the sum of squared errors
% between the frequency points

% based on model that s = m.(r.f) + c {+a.((r.f).^2)} (not implemented)

% alignment could be weighted to the intensity of the peaks (ie multiply
% the error by the data)

% INPUTS:
% rL,sL,rR,sR: reference and subject (ie subject) spectra, containing
%   fields .y and .f: left and right (as applicable).  Align only non-empty
%   spectra
% params: [c,m,a] offset (0 order), scaling (1st order) and 2nd order
%   (not implemented) parameters. m and a parameters are optional

% Change log

% 03/Mar/08 v0.1    Initial version

c = params(1);
if length(params) == 2
    m = params(2) / 1e6;    % to improve sensitivity of fminsearch
else
    m = 1;
end

% number of peaks using to align
if isempty(rL)
    nL = 0;
else
    nL = length(rL.f);
end
if isempty(rR)
    nR = 0;
else
    nR = length(rR.f);
end

% apply alignment
if nL
    sL.f = m*sL.f + c;
end
if nR
    sR.f = m*sR.f + c;
end

% sum of squares
resid = 0;
if nL
    resid = resid + sum( (sL.f-rL.f).^2 );
end
if nR
    resid = resid + sum( (sR.f-rR.f).^2 );
end


function [SDout, SPout] = SegStitch_1_1(SD, SP, P, marginToCut)
%stitches multiple spectra together, assuming A and B are the same for all
%spectra, using spline interpolation to make the data points monotonic
%assumptions:
%noise has been removed from the spectra (ie set to zero)
%limitations:
%the frequency points must be evenly spaced
%the spacing in the first range will determine the spacing of the stitched spectrum

% INPUTS:
% SD: array of spectra data structs, one for each SIM window
% SP: array of spectra peak structs, one for each SIM window
% P: array of parameters structs, one for each SIM window
% marginToCut: edge effects in each SIM window to remove during stitching

% OUTPUTS:
% SDout: single spectra data struct (stitched)
% SPout: single spectra peak struct (stitched)

%version history:
%9/5/07, v0.2    redundant modes removed (leaving spline only)
%                removed references to specdf
%18/6/07         dead region now dependent upon window mzStart
%19/6/07         included test if no overlapping regions left, in which case peaks are concat to save time
%                added check to remove any data points beyond the spectral mzStart or mzEnd
%06/02/08, v0.3  replaced spec.peaksn with spec.peaksFlag and spec.pNoise
%                only cut spectra at minima below the data threshold
%22/02/08        minor error-causing bug fix
%26/02/08, v1.0  Removed params, numSpec and mode parameters, tidied code
% 13/Mar/08         Added output of by how much target cut point is moved
% 14/Apr/08         Minor change to remove old specParams variable line 78
% 14/Apr/08 v1.1    Added edge effect to remove from start of first SIM
%                   window and end of last SIM window

numSpec = length(P);

% catch case only 1 SIM window
if numSpec == 1
    SDout = SD;
    SPout = SP;
    return
end

% determine a common freq. for each spectrum, will form the freq points of the new stitched spectrum
% the deltaf should be taken from a CALIBRATED (ie freq. points have not been adjusted) spectrum
deltaf = 0;
for si = 1:numSpec
    if (P(si).ecald || P(si).icald)
        deltaf = SD(si).f(2) - SD(si).f(1);
        break;
    end
end
if deltaf == 0
    error('No calibrated SIM windows found');
end

%cut the "edge effect" from spectra
%select edge effect based on midpoint of overlap
sw_l = zeros(numSpec,1);
sw_h = zeros(numSpec,1);
sw_l(1) = marginToCut.MZMIN;
for si=2:numSpec
    mp = (P(si).mzStart + P(si-1).mzEnd)/2;
    for i=1:length(marginToCut.BOUND)-1
        if (marginToCut.BOUND(i) < mp && mp <= marginToCut.BOUND(i+1))
            sw_l(si) = marginToCut.START(i);
            sw_h(si-1) = marginToCut.END(i);
        end
    end
end
sw_h(end) = marginToCut.MZMAX;
for si=1:numSpec
    fprintf('Segment %d (%d-%dm/z): %.2f (lower edge), %.2f (higher edge)\n',...
        si,P(si).mzStart,P(si).mzEnd,sw_l(si),sw_h(si));
    %if after removing the edge effects the windows still overlap, unpreferentially increase the edge effects until no overlap exists
    %this is to avoid averaging spectral values
    if si<numSpec
        olap = P(si).mzEnd - P(si+1).mzStart - sw_h(si) - sw_l(si+1);
        if olap > 0
            fprintf('INCREASING AMOUNT OF REMOVED EDGE EFFECT to zero overlap:\n');
            sw_h(si) = sw_h(si) + olap /2;
            sw_l(si+1) = sw_l(si) + olap /2;
            fprintf('Segment %d (%d-%dm/z): %.2f (lower edge), %.2f (higher edge)\n',si,P(si).mzStart,P(si).mzEnd,sw_l(si),sw_h(si));
        elseif olap < 0
            error('Too much edge effect: missing data points');
        end
    end
    if P(si).dNthresh > 0
        %only cut at a minima in the data below the threshold (if noise
        %non-zero)
        zIdx = maxima(-1*SD(si).y,length(SD(si).y)) + 1;
        zIdx = zIdx(SD(si).y(zIdx) < P(si).dNthresh);
    elseif length(find(SD(si).y==0)) > 0
        % only cut at a zero (if data contains zero values)
        zIdx = find(SD(si).y==0);
    else
        % only cut at a minimum (if noise not set and no zeros)
        zIdx = maxima(-1*SD(si).y,length(SD(si).y)) + 1;
    end
    %higher frequency end
    fHiCut = mz2f(P(si).mzStart + sw_l(si),P(si));
    idx = find(SD(si).f <= fHiCut);
    iH = max(intersect(idx, zIdx));
    % check haven't moved way beyond target cut point
    fprintf('Moved low cut point by %.2f m/z\n',...
        abs(f2mz(SD(si).f(max(idx)),P(si))-f2mz(SD(si).f(max(iH)),P(si))));
    %lower frequency end
    fLoCut = mz2f(P(si).mzEnd - sw_h(si),P(si));
    idx = find(SD(si).f >= fLoCut);
    iL = min(intersect(idx, zIdx));
    % check haven't moved way beyond target cut point
    fprintf('Moved high cut point by %.2f m/z\n',...
        abs(f2mz(SD(si).f(min(idx)),P(si))-f2mz(SD(si).f(max(iL)),P(si))));
    %cut
    sFields = fieldnames(SD);
    for i=1:length(sFields)
        SD(si).(sFields{i}) = SD(si).(sFields{i})(iL:iH);
    end

    %peaks - high freq
    iH = max(find(SP(si).f <= fHiCut));
    %lower frequency end
    iL = min(find(SP(si).f >= fLoCut));
    %cut
    sFields = fieldnames(SP);
    for i=1:length(sFields)
        SP(si).(sFields{i}) = SP(si).(sFields{i})(iL:iH);
    end
end

% trim spectra to mzStart and mzEnd if necessary
for si=1:numSpec
    % spectral data
    mz = f2mz(SD(si).f,P(si));
    idx = find(P(si).mzEnd < mz | mz < P(si).mzStart);
    if ~isempty(idx)
        fprintf('Trimming spec %d (%.2f-%.2f): removing spec data from %.6f-%.6f\n',...
            si,P(si).mzStart,P(si).mzEnd,min(mz(idx)),max(mz(idx)));
        sFields = fieldnames(SD);
        for i=1:length(sFields)
            SD(si).(sFields{i})(idx) = [];
        end
    end
    % spectral peaks
    mz = f2mz(SP(si).f,P(si));
    idx = find(P(si).mzEnd < mz | mz < P(si).mzStart);
    if ~isempty(idx)
        fprintf('Trimming spec %d (%.2f-%.2f): removing peaks from %.6f-%.6f\n',...
            si,P(si).mzStart,P(si).mzEnd,min(mz(idx)),max(mz(idx)));
        sFields = fieldnames(SP);
        for i=1:length(sFields)
            SP(si).(sFields{i})(idx) = [];
        end
    end
end

%next determine what Sf, the f values for the final spectrum, will be
%assuming the frequency values are already cut
%find the maximum and minimum values of Sf
SfMin = min(SD(1).f);
SfMax = max(SD(1).f);
for si = 1:numSpec
    fMax(si) = max(SD(si).f);
    fMin(si) = min(SD(si).f);
    SfMax = max([SfMax fMax(si)]);
    SfMin = min([SfMin fMin(si)]);
end

%create Sf
Sf = SfMax:-deltaf:SfMin;
%create max to min then flip to avoid a /0 error
%also makes sense to construct based on highest freq. (1st range).
Sf = fliplr(Sf);

% 5/12/06 - bug fix... if the first segment has been aligned, need to remove the offset here such that the
% first segment (because stitching is referenced from the 1st seg) aligns with original data points for
% KCe interpolation to work.
m = Sf(1) / deltaf;
segShift = (round(m)-m)*deltaf;
Sf = Sf + segShift;

%if no regions are now overlapping, simply concatenate the remaining
%lengths of windows (spectrum may not have monotonous frequency values for
%peak-picking)
if issorted([SD(end:-1:1).f])
    fprintf('No overlap remaining: concatenating windows...');
    % spectral data
    sFields = fieldnames(SD);
    for i=1:length(sFields)
        SDout.(sFields{i}) = [SD(end:-1:1).(sFields{i})];
    end
    % spectral peaks
    sFields = fieldnames(SP);
    for i=1:length(sFields)
        SPout.(sFields{i}) = [SP(end:-1:1).(sFields{i})];
    end
    %check consistency
    if ~issorted(SDout.f), error('Inconsistency in stitching'); end
    if ~issorted(SPout.f), error('Inconsistency in stitching'); end
    fprintf('done\n');
else
    error('Code not capable of stitching overlapping windows');
end
