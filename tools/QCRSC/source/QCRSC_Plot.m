function QCRSC_Plot(html_indir, html_outfile, html_outdir, peaks)

mkdir(html_outdir);

config_file = [html_indir, filesep, 'configx.txt'];

config = fileread(config_file);
config = strrep(config, 'data.csv', [html_indir, filesep, 'data.csv']);
config = strrep(config, 'peaks.csv', [html_indir, filesep, 'peaks.csv']);
config = strrep(config, '%', '%%');

config_file_temp = tempname;
fileID = fopen(config_file_temp, 'w'); %Make sure paths are pointing to correct path for plotting tool
fprintf(fileID, config);
fclose(fileID);

poly_file = [html_indir, filesep, 'Filt_output_data.csv'];
BeforeAfter_file = [html_indir, filesep, 'QCRSC_output_data.csv'];

fid = fopen(html_outfile, 'wt');
fprintf(fid,'<?xml version="1.0" encoding="utf-8" ?>');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
fprintf(fid,'<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">');
fprintf(fid,'<head>');
%fprintf(fid,'<title>Sum Transient Data - Output</title>');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />');
fprintf(fid,'<link rel="stylesheet" href="/static/style/base.css" type="text/css" />');
fprintf(fid,'</head>');
fprintf(fid,'<body>');
fprintf(fid,'<div class="document">');
fprintf(fid,'<div class="donemessagelarge">');
fprintf(fid,'<div style="padding: 3px"><h2><b>Figures - QCRSC</b></h2></div>');
fprintf(fid,'<hr></hr>');


if exist(config_file, 'file') && exist(poly_file, 'file')
    fprintf(fid,'<div style="padding: 3px"><h3><b>Peak Outlier Detection</b></h2></div>');
    [fns_poly] = PlotPolyFilterData(config_file_temp, peaks, html_outdir);
    for i=1:length(fns_poly)
        if ~isempty(fns_poly{i})
            fn = strrep(fns_poly{i},[html_outdir,filesep],'');
            fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn,'" target="_blank">',fn,'</a></b></div>']);
        end
    end
end

if exist(config_file, 'file') && exist(BeforeAfter_file, 'file')
    [fns_QCRSC] = PlotQCRSCData(config_file_temp, peaks, html_outdir);
    [fns_BeforeAfter] = PlotBeforeAfter(config_file_temp, peaks, html_outdir);
    fprintf(fid,'<div style="padding: 3px"><h3><b>Peak within batch ..</b></h2></div>'); 
    for i=1:length(fns_QCRSC)
        if ~isempty(fns_QCRSC{i})
            fn = strrep(fns_QCRSC{i},[html_outdir,filesep],'');
            fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn,'" target="_blank">',fn,'</a></b></div>']);
        end
    end
    fprintf(fid,'<div style="padding: 3px"><h3><b>Peak accross batches (Before and after batch correction)</b></h2></div>'); 
    for i=1:length(fns_BeforeAfter)
        if ~isempty(fns_BeforeAfter{i})
            fn = strrep(fns_BeforeAfter{i},[html_outdir,filesep],'');
            fprintf(fid, ['<div style="padding: 3px"><b><a href="',fn,'" target="_blank">',fn,'</a></b></div>']);
        end
    end
end

fprintf(fid,'<hr></hr>');
fprintf(fid,'</div>');
fprintf(fid,'</div>');
fprintf(fid,'</body></html>');
fclose(fid);






