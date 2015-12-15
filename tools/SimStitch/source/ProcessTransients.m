function ProcessTransients(fileList, html_indir, html_outfile, html_outdir)
%
%
% ProcessTransients_BB(fileList,html_indir, html_outfile, html_outdir)
%
% The 2nd step in SimStitch processing of DIMS FTICR data.
%
% inputs:
%	fileList:	the full path of an xml file containing information regarding location of the RAW files etc and some communal variables
%	html_indir: 	the full path to the directory containing the files referenced in 'html_infile'. This should replace the averaged_transients_directory variable in fileList.  
%	html_outfile:	the full path for a new html file that will reference the files output by ProcessTransients for use by Galaxy.
%	html_outdir:	the full path for a directory in which to store the files created by ProcessTransients that will be referenced in html_outfile by Galaxy
%
% outputs:
%	this function does not return any variables for Matlab because it is designed to create output files as expected by Galaxy. 
%
%
%
% This file has been adapted for Galaxy from many previous incarnations by R.L.Davidson
% For more information contact Mark Viant m.viant@bham.ac.uk 
%
%
% 
% 23/04/2015
% JE: Cleaned up code; added support for ... 

%set input parameters to fit existing legacy code
QUIET = 1;
ZFILLS = 1;

%extract the fileList information from the xml file
fileListpath = fileList;
fileList = ImportFileListXML(fileList);

%make directory for processed transients
fileList.ATDir = [html_indir,filesep];
fileList.FSDir = [html_outdir,filesep]; %this converts the path given by Galaxy for the output files, adding a forward slash (may need changed for other OS)

dFold = fileList.FSDir;
try mkdir(dFold);
catch 
    error(['Cannot make processed transient directory ',dFold]); 
end

filename_html = {}; %this is for storing Galaxy output

%get parameter information for each spectrum to be stitched
for fi=1:fileList.nDataFiles
    fprintf('Processing transient files, stem: %s\n',fileList.Samples(1,fi).ID);
    FS = [];
    FSParams = [];
    isSingle = 0;   %don't change
    isStitch = 0;   %don't change
    listFilterIndices = [];
    f = 1;  %filter index
    while 1
        %load from matlab files the averaged transient(s)
        if ~isStitch && ~isSingle
            try
                load(fileList.Samples(1,fi).ID); %single scan
                isSingle = 1;
            catch
                try
                    load([fileList.ATDir,fileList.Samples(1,fi).ID,'_seg01.mat']);
                    disp('_seg01.mat');
                    isStitch = 1;
                catch
                    error(['file ',fileList.Samples(1,fi).ID,' or ',fileList.ATDir,fileList.Samples(1,fi).ID,'_seg01.mat not found']);
                end
            end
        else
            if isStitch
                try
                    load([fileList.ATDir,fileList.Samples(1,fi).ID,'_seg',num2str(f,'%.2d'),'.mat']);
                    fprintf('_seg%.2d.mat\n',f);
                catch
                    %no more files in the stitch series
                    break;
                end
            else
                %single scan and done
                break;
            end
        end

        if exist('specParamsT'), FSParams = specParamsT; clear specParamsT; end

        %correct for missing variables
        if ~isfield(FSParams,'filterIndex')
            warning('filterIndex missing from averaged transients; acquiring from raw file');
            FSParams.filterIndex = f;
        end

        listFilterIndices = unique([listFilterIndices FSParams.filterIndex]);

        FSParams_f = FSParams([FSParams.filterIndex]==listFilterIndices(f));

        %zero-fills
        FSParams_f.zFills = ZFILLS;

        %sum the TIC over scans for this range
        FSParams_f.TIC = sum(FSParams_f.TIC);

        %maximise potential precision (16-bit limits to m/z of 0.2ppm @ 70m/z)
        transient = double(transient);

        % Get the frequency spectrum (replaces ProcessTrans function from old DIMS pipeline code)
        % a. remove DC offset
        transient = transient - mean(transient);
        % b. perform apodisation
        transient = Apodize(transient,'hanning',QUIET);
        % c.convert to freq. domain
        n = length(transient);
        % TO DO: Add support for bruker machine
        f1 = mz2f(FSParams_f.mzEnd,FSParams_f); %lower freq.
        f2 = mz2f(FSParams_f.mzStart,FSParams_f); %upper freq.
        % Zero-fill, fourier-transform and calculate magnitude spectrum
        % (replaces Transform function)
        % TO DO: Support for bruker machine; upper lower frequency
        % limit
        n = (2^ZFILLS) * n; % Double number of data points for zero filling
        FS(f).data = fft(transient, n);     %fft'd data
        FS(f).data = FS(f).data(1:n/2+1);     %first half of spectrum required
        fs = 2*FSParams_f.BW; %Sampling frequency
        df = fs/n;  %freq step
        FS(f).f = (0:n/2)*df;    %frequency scale
        %chop frequency
        firstIdx = find(FS(f).f>f1, 1 );
        lastIdx = find(FS(f).f<f2, 1, 'last' );
        FS(f).f = FS(f).f(firstIdx:lastIdx); %Truncate frequency scale
        FS(f).data = FS(f).data(firstIdx:lastIdx); %Truncate frequency spectrum
        FS(f).data = abs(FS(f).data); %magnitude mode      
        % d.baseline correction (optional; code was doing nothing so it was
        % removed)

        %save the params (filter)
        FSParams_all(f) = FSParams_f;

        %next filter
        f = f+1;

    end

    %rename parameters variable
    FSParams = FSParams_all;
    clear FSParams_all FSParams_f transient;

    %the number of SIM windows
    numWindows = length(FSParams);

    %check data and freq's and params same length
    for si=1:numWindows
        if length(FS(si).data)~=length(FS(si).f), error('FS(i).data and FS(i).f size mismatch'); end
    end
    if length(FS)~=length(FSParams), error('FS and FSParams size mismatch'); end        

    %save spectrum
    fileName = [fileList.FSDir,'FS','_',fileList.Samples(1,fi).ID];
    fprintf('Saving frequency spectrum to %s.mat...',fileName);
    save(fileName,'FS','FSParams');
    fprintf('done\n');
    
    %prepare to save to html_outfile for Galaxy
    filename_html{end+1} = ['FS','_',fileList.Samples(1,fi).ID,'.mat'];

end
%create html output for Galaxy
fid = fopen(html_outfile, 'wt');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
fprintf(fid,'<html><head>');
fprintf(fid,'<title>Processed Transients - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">');
fprintf(fid,'<link href="/static/style/base.css?v=1415642706" media="screen" rel="stylesheet" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'   <div style="padding: 3px"><h2><b>Frequency Spectra</b></h2></div>');
fprintf(fid,'<hr></hr>');
for i=1:length(filename_html)
  	html_code = ['<div style="padding: 3px"><b><a href="',filename_html{i},'">',filename_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);
