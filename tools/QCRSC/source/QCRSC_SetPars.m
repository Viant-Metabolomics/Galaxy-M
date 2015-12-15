function QCRSC_SetPars(dso_xml, html_outfile, html_outdir, toolname, varargin)

mkdir(html_outdir);

[pathstr,~,~] = fileparts(mfilename('fullpath'));

config = fileread([pathstr, filesep, 'configx.txt']);

config = strrep(config, '$output_files_path', html_outdir);
if strcmp(toolname, 'poly')
	config = strrep(config, '$runPolyFilter', 'yes');
	config = strrep(config, '$runQCRSC', 'no');
	config = strrep(config, '$runPQF', 'no');
	config = strrep(config, '$peakTransform', varargin{1});
	config = strrep(config, '$polyFunc', varargin{2});
	config = strrep(config, '$polyCI', varargin{3});
	config = strrep(config, '$polyAction', varargin{4});
	config = strrep(config, '$QCRSCoperator', 'subtract'); %default value - not used
	config = strrep(config, '$QCRSCmaxWindow', '-1'); %default value - not used
	config = strrep(config, '$QCRSCsearchRangeMin', '0'); %default - not used
	config = strrep(config, '$QCRSCsearchRangeSteps', '0.5'); %default - not used
	config = strrep(config, '$QCRSCsearchRangeMax', '7'); %default - not used
	config = strrep(config, '$KsiTol', '4'); %default value - not used
	config = strrep(config, '$kill', 'yes'); %default value - not used
	config = strrep(config, '$pdiff', '1e-4'); %default value - not used
	config = strrep(config, '$QCdist', '3'); %default value - not used
	config = strrep(config, '$maxRSD', '25'); %default value - not used
	config = strrep(config, '$minD_RATIO', '1'); %default value - not used
elseif strcmp(toolname, 'BC')
	config = strrep(config, '$runPolyFilter', 'no');
	config = strrep(config, '$runQCRSC', 'yes');
	config = strrep(config, '$runPQF', varargin{1});
	config = strrep(config, '$peakTransform', varargin{2});
	config = strrep(config, '$polyFunc', 'poly3'); %default value - not used
	config = strrep(config, '$polyCI', '0.99'); %default value - not used
	config = strrep(config, '$polyAction', 'ignore'); %default value - not used
	config = strrep(config, '$QCRSCoperator', varargin{3});
	config = strrep(config, '$QCRSCmaxWindow', varargin{4});
	config = strrep(config, '$QCRSCsearchRangeMin', varargin{5});
	config = strrep(config, '$QCRSCsearchRangeSteps', varargin{6});
	config = strrep(config, '$QCRSCsearchRangeMax', varargin{7});
	config = strrep(config, '$KsiTol', varargin{8});
	config = strrep(config, '$kill', varargin{9});
	config = strrep(config, '$pdiff', varargin{10});
	config = strrep(config, '$QCdist', varargin{11});
	config = strrep(config, '$maxRSD', varargin{12});
	config = strrep(config, '$minD_RATIO', varargin{13});
else
	fprintf('Check number of variables');
end

config = strrep(config, '%', '%%');
config

config_file = [html_outdir, filesep, 'configx.txt'];
fileID = fopen(config_file, 'w');
fprintf(fileID, config);
fclose(fileID);

csv_data = [html_outdir, filesep, 'data.csv'];
csv_peaks = [html_outdir, filesep, 'peaks.csv'];

%copyfile([pathstr, filesep, 'data.csv'], [html_outdir, filesep]); %DEBUG
%copyfile([pathstr, filesep, 'peaks.csv'], [html_outdir, filesep]); %DEBUG

dso2csv(dso_xml, csv_data, csv_peaks)
QCRSC_X3(config_file)

fid = fopen(html_outfile, 'wt');
fprintf(fid,'<?xml version="1.0" encoding="utf-8" ?>');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
fprintf(fid,'<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">');
fprintf(fid,'<head>');
%fprintf(fid,'<title></title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />');
fprintf(fid,'<link rel="stylesheet" href="/static/style/base.css" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="document">');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'<div style="padding: 3px"><h2><b>Outputs - QCRSC</b></h2></div>');
fprintf(fid,'<hr></hr>');


fprintf(fid,'<div style="padding: 3px"><h3><b>Configuration and Input Files</b></h2></div>');
fprintf(fid, '<div style="padding: 3px"><b><a href="configx.txt" target="_blank">configx.txt</a></b></div>');
fprintf(fid, '<div style="padding: 3px"><b><a href="peaks.csv" target="_blank">peaks.csv</a></b></div>');
fprintf(fid, '<div style="padding: 3px"><b><a href="data.csv" target="_blank">data.csv</a></b></div>');

fn_peaks_Filt = 'Filt_output_peaks.csv';
fn_data_Filt = 'Filt_output_data.csv';
fn_peaks_QCRSC = 'QCRSC_output_peaks.csv';
fn_data_QCRSC = 'QCRSC_output_data.csv';
fn_peaks_Cleaned = 'Cleaned_output_peaks.csv';
fn_data_Cleaned = 'Cleaned_output_data.csv';

varargin

if strcmp(toolname, 'poly')
        
    csv2dso(strrep(csv_data, 'data.csv', fn_data_Filt), strrep(csv_peaks, 'peaks.csv', fn_peaks_Filt), varargin{5})
    
    fprintf(fid,'<div style="padding: 3px"><h3><b>Peak Outlier Detection</b></h2></div>');
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_data_Filt,'" target="_blank">',fn_data_Filt,'</a></b></div>']);
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_peaks_Filt,'" target="_blank">',fn_peaks_Filt,'</a></b></div>']);

elseif strcmp(toolname, 'BC') && strcmp(varargin{14}, 'None')
    
    csv2dso(strrep(csv_data, 'data.csv', fn_data_QCRSC), strrep(csv_peaks, 'peaks.csv', fn_peaks_QCRSC), varargin{14})
    
    fprintf(fid,'<div style="padding: 3px"><h3><b>Batch Correction</b></h2></div>');
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_data_QCRSC,'" target="_blank">',fn_data_QCRSC,'</a></b></div>']);
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_peaks_QCRSC,'" target="_blank">',fn_peaks_QCRSC,'</a></b></div>']);

else
    csv2dso(strrep(csv_data, 'data.csv', fn_data_QCRSC), strrep(csv_peaks, 'peaks.csv', fn_peaks_QCRSC), varargin{14});
    csv2dso(strrep(csv_data, 'data.csv', fn_data_Cleaned), strrep(csv_peaks, 'peaks.csv', fn_peaks_Cleaned), varargin{15});

    fprintf(fid,'<div style="padding: 3px"><h3><b>Batch Correction</b></h2></div>');
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_peaks_QCRSC,'" target="_blank">',fn_peaks_QCRSC,'</a></b></div>']);
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_data_QCRSC,'" target="_blank">',fn_data_QCRSC,'</a></b></div>']);
    fns = dir(html_outdir);
    fns =  {fns.name};
    for i=1:length(fns)
        if ~isempty(strfind(fns{i}, '_Batch_'))
            fprintf(fid, ['<div style="padding: 3px"><b><a href="',fns{i},'" target="_blank">',fns{i},'</a></b></div>']);
        end
    end
    fprintf(fid,'<div style="padding: 3px"><h3><b>Spectral Cleaning</b></h2></div>');
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_peaks_Cleaned,'" target="_blank">',fn_peaks_Cleaned,'</a></b></div>']);
    fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn_data_Cleaned,'" target="_blank">',fn_data_Cleaned,'</a></b></div>']);
    fprintf(fid, '<div style="padding: 3px"><b><a href="output_Removed.csv" target="_blank">output_Removed.csv</a></b></div>');
end

fprintf(fid,'<hr></hr>');
fprintf(fid,'</div>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);

fileID = fopen(config_file, 'w'); %Make sure paths are pointing to correct path for plotting tool
fprintf(fileID, strrep(config, [html_outdir, filesep], ''));
fclose(fileID);

