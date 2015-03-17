function ProcessTransients_BB(fileList, html_infile, html_indir, html_outfile, html_outdir)
%
%
% ProcessTransients_BB(fileList,html_infile, html_indir, html_outfile, html_outdir)
%
% The 2nd step in SimStitch processing of DIMS FTICR data.
%
% inputs:
%	fileList:t	the full path of an xml file containing information regarding location of the RAW files etc and some communal variables
%	html_infile:	the full path of an html file used by Galaxy to allow the user to view and download multiple files output by any tool. This HTML file 'could' be parsed to discover filenames but fileList currently fulfils that role; a legacy from previous SimStitch incarnations that avoids HTML parsing here.
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

%set input parameters to fit existing legacy code
SET = fileList; %this needs CHANGED now that an individual xml file is being passed



QUIET = 1;
ZFILLS = 1;

%extract the fileList information from the xml file
fileList = import_filelist_xml(SET);

fileList.scans = {};

%make directory for processed transients
fileList.avTransDir = [html_indir,filesep]; %this updates the location of the files output by SumTransients to reflect the location utilised by Galaxy
fileList.specDir = [html_outdir,filesep]; %this converts the path given by Galaxy for the output files, adding a forward slash (may need changed for other OS)

dFold = fileList.specDir;
try mkdir(dFold);
catch 
    error(['Cannot make processed transient directory ',dFold]); 
end

rawCount = length(fileList.rawFull);

filename_html = {}; %this is for storing Galaxy output

%get parameter information for each spectrum to be stitched
for fi=1:rawCount
    fprintf('Processing transient files, stem: %s\n',fileList.rawShort{fi});
    spec = [];
    specParams = [];
    isSingle = 0;   %don't change
    isStitch = 0;   %don't change
    listFilterIndices = [];
    f = 1;  %filter index
    while 1
        %load from matlab files the averaged transient(s)
        if ~isStitch & ~isSingle
            try
                load(fileList.rawShort{fi}); %single scan
                isSingle = 1;
            catch
                try
                    load([fileList.avTransDir,fileList.rawShort{fi},'_seg01.mat']);
                    disp('_seg01.mat');
                    isStitch = 1;
                catch
                    error(['file ',fileList.rawShort{fi},' or ',fileList.avTransDir,fileList.rawShort{fi},'_seg01.mat not found']);
                end
            end
        else
            if isStitch
                try
                    load([fileList.avTransDir,fileList.rawShort{fi},'_seg',num2str(f,'%.2d'),'.mat']);
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

        if exist('specParamsT'), specParams = specParamsT; clear specParamsT; end

        %correct for missing variables
        if ~isfield(specParams,'filterIndex')
            warning('filterIndex missing from averaged transients; acquiring from raw file');
            specParams.filterIndex = f;
        end

        listFilterIndices = unique([listFilterIndices specParams.filterIndex]);
%             warning('modified filter handling for sliding seg expt');
%             listFilterIndices = [listFilterIndices specParams.filterIndex];

        specParams_f = specParams(find([specParams.filterIndex]==listFilterIndices(f)));

        %zero-fills
        specParams_f.zFills = ZFILLS;

        %sum the TIC over scans for this range
        specParams_f.TIC = sum(specParams_f.TIC);

% % USE ONLY PART OF THE TRANSIENT  - INVESTIGATION 14/May/07
%             Td = 10e-3;      %window time duration (s)
%             Ts = 0;      %window start time (s)
%             transient = transient(floor(Ts / (specParams.T/length(transient))+1):floor((Ts + Td) / (specParams.T/length(transient))));
%             specParams_f.T = Td;
% % END

% TEST 01/June/07 - get certain ranges only
%             if 70>specParams_f.mzStart | specParams_f.mzStart>80, f=f+1; continue;
%             else specParams_f.mzStart = 50; specParams_f.mzEnd = 150; fprintf('.'); end
% END TEST

        %maximise potential precision (16-bit limits to m/z of 0.2ppm @ 70m/z)
        transient = double(transient);

        %get the spectrum
        [spec(f).data, spec(f).f] = ProcessTrans(transient,specParams_f, QUIET);

        %save the params (filter)
        specParams_all(f) = specParams_f;

        %next filter
        f = f+1;

    end

    %rename parameters variable
    specParams = specParams_all;
    clear specParams_all specParams_f transient;

    %the number of spectra
    numSpectra = length(specParams);

    %check data and freq's and params same length
    for si=1:numSpectra
        if length(spec(si).data)~=length(spec(si).f), error('spec(i).data and spec(i).f size mismatch'); end
    end
    if length(spec)~=length(specParams), error('spec and specParams size mismatch'); end        

    %save spectrum
    fileName = [fileList.specDir,'spec',fileList.setUID,'_',fileList.rawShort{fi}];
    fprintf('Saving overlapping spectra to %s.mat...',fileName);
    save(fileName,'spec','specParams');
    fprintf('done\n');
    
    %prepare to save to html_outfile for Galaxy
    filename_html{end+1} = ['spec',fileList.setUID,'_',fileList.rawShort{fi},'.mat'];

end
%create html output for Galaxy
fid = fopen(html_outfile, 'wt');
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
for i=1:length(filename_html)
  	html_code = ['<div style="padding: 3px"><b><a href="',filename_html{i},'">',filename_html{i},'</a></b></div>'];
	fprintf(fid, html_code);
end
fprintf(fid,'<hr></hr>');
fprintf(fid,'<p>');
fprintf(fid,'Note: ---');
fprintf(fid,'</p>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);
return


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



function [data,f] = ProcessTrans(transient,params,QUIET)

% MODIFICATIONS
% 
% 13/Mar/08     Updated mz2f to new function definition



PLOT_ON = 0;

%remove DC offset
if ~QUIET, disp('Removing DC offset...'); end
transient = transient - mean(transient);
if ~QUIET, disp('  ...done'); end
if PLOT_ON, plot(transient,'k'); pause; end   %f in MHz

%perform apodisation
transient = Apodize(transient,'hanning',QUIET);
if PLOT_ON, plot(transient,'k'); pause; end   %f in MHz

%convert to freq. domain
n = length(transient);
f1 = mz2f(params.mzEnd,params); %lower freq.
f2 = mz2f(params.mzStart,params); %upper freq.
[data,f] = Transform(transient,params.BW,'fft','magnitude',params.zFills,f1,f2,n,QUIET);
if PLOT_ON, plot(f./1e6,data); hold all; axis('auto'); pause; end   %f in MHz

%baseline correction
data = BaseLine(data,'off', QUIET);
return



function vOutput = Apodize(vInput,method,QUIET)
% This function takes a vector input and outputs the same vector,
% multiplied by the apodization window function



%Determine the number of input data points
n = length(vInput);

%Initialize the vector
vFunc = linspace(0,n-1,n);

switch method
    case {0,'none'}
        if ~QUIET, disp('No apodisation'); end
        vOutput = vInput;
        return
    case {1,'hanning'}
        %Calculate the hanning funtion
        if ~QUIET, disp('Applying Hanning...'); end
        vFunc = .5*(1-cos(2*pi*vFunc/(n-1)));
    case {1,'hamming'}
        %Calculate the hamming funtion
        if ~QUIET, disp('Applying Hamming...'); end
        vFunc = 0.54 - 0.46*cos(2*pi*vFunc/(n-1));
    case {3,'blackman-harris'}
        if ~QUIET, disp('Applying Blackmann-Harris...'); end
        vFunc = 0.42323 - 0.49755*cos(2*pi*vFunc/(n-1)) + 0.07922*cos(4*pi*vFunc/(n-1));
    otherwise
        disp('ERROR: Choice not recognised! Default: no apodisation selected.');
        vOutput = vInput;
        return
end

%Output the result
vOutput = vInput.*vFunc;
if ~QUIET, disp('  ...done'); end
return



function [fdata,f] = Transform(tdata, BW, transform, mode, zeroFill, f1, f2, n, QUIET)
%Perform frequency domain transformation with options



switch transform
    case {0, 'fft'}
        if isempty(zeroFill)
            disp('WARNING: No zero-fill specified, Default: 1 zero-fill selected');
            zeroFill = 1;
        end
        if ~QUIET, disp(['Applying FFT with ',num2str(zeroFill),' zero fill(s), BW=',num2str(BW),'...']); end
        tic;
        n = (2^zeroFill) * length(tdata);
        fdata = fft(tdata, n);     %fft'd data
        fdata = fdata(1:n/2+1);     %first half of spectrum required
        fs = 2*BW;
        df = fs/n;  %freq step
        f = (0:n/2)*df;    %frequency scale
        %mzEps = 0;      %add m/z 1 margin overlap
        %f1 = f1 - mz2f(mzEps,A,B);
        %f2 = f2 + mz2f(mzEps,A,B);
        %chop frequency
        firstIdx = min(find(f>f1));
        lastIdx = max(find(f<f2));
        f = f(firstIdx:lastIdx);
        fdata = fdata(firstIdx:lastIdx);
    case {1, 'chirp'}
        disp(['Applying CZT with BW=',num2str(BW),'...']);
        tic;
        fs = BW*2;
        w = exp(-j*2*pi*(f2-f1)/(n*fs));
        a = exp(j*2*pi*f1/fs);
        fdata = czt(tdata,n,w,a);
        f = ((0:length(fdata)-1)*(f2-f1)/length(fdata)) + f1;
    otherwise
        disp('ERROR: Transform choice not recognised!');
        return
end

if ~QUIET, disp([' ...done (elapsed time: ',num2str(toc),'s)']); end

switch mode
    case {0,'magnitude'}
        if ~QUIET, disp('Calculating magnitude-mode...'); end
        fdata = abs(fdata); %magnitude mode
    otherwise
        disp('ERROR: Mode choice not recognised!');
        return
end

if ~QUIET, disp('  ...done'); end
return



function outData = BaseLine (inData, mode, QUIET)
%Baseline correct


switch mode
    case {0,'off'}
        if ~QUIET, disp('No baseline correction'); end
        outData = inData;
        return
    %case {1, 'on'}
        %disp('Applying baseline correction...');
    otherwise
        disp('ERROR: Mode choice not recognised!');
        outData = inData;
        return
end

if ~QUIET, disp('  ...done'); end
return


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
return



