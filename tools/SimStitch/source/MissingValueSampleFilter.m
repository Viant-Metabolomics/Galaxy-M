function [num_missing, stats, flag_samples] = MissingValueSampleFilter(dso_xml, sd_flag, tails, conf, del, figurefile, outfile)

%
% num_missing = MissingVals(dso, sd_flag, tails, conf)
%
%   inputs:
%           dso: a dataset object (pls toolbox)
%           sd_flag: 1 or 0 to use standard deviation or standard error
%           respectively
%           tails: 1 for 1-sided, 0 for 2-sided
%           conf: the percentage confidence you want e.g. 95
%
%   outputs:    
%           num_missing: a vector where each entry represents a sample from
%           the dso input. 
%           stats: a struct with the following variables
%               .mu         - the mean of missing values
%               .upper      - the upper estimate
%               .lower      - the lower estimate, only if 2-sided
%               .tails      - records 1-sided or 2-sided
%               .method     - standard deviation or standard error
%               .conf       - the confidence level 
%               .df         - the degrees of freedom (number of samples-1)
%               .z          - number of standard deviations/errors used to
%                           calculate the limits 
%                           e.g. 95% confidence with 2-sides uses
%                           1.96 standard deviations when there are
%                           infinite samples but uses 2.04 when there are 
%                           only 30 samples.
%                   
%
%   displays: 
%           plots a bar chart of the missing values.
%
%
%   r.l.davidson 19/03/2013
%
%   Version 1.2

[dso, name, source] = autoimport(dso_xml, 'xml');

if isa(sd_flag, 'char')
	sd_flag = str2num(sd_flag);
end

if isa(tails, 'char')
	tails = str2num(tails);
end

if isa(conf, 'char')
	conf = str2num(conf);
end

if isa(del, 'char')
	del = str2num(del);
end


sam_include = dso.include{1};
if isempty(sam_include)
    disp('error: no samples included')
    num_missing = [];
    return
end

peak_include = dso.include{2};
data_orig = dso.data(sam_include,peak_include);

data_bin = data_orig;		%make a new copy for �binary� values

data_bin(find(data_orig>0))=1;		%set all values greater than 0 to �1�

data_bin = data_bin-1;		%now all zeros = -1 and all non-zeros = 0;

num_missing = abs(sum(data_bin'));		%then we sum up the missing values
n = length(num_missing);  
m = size(data_bin,2);

perc_missing = num_missing./m*100;   %convert to percentage missing

sd = std(perc_missing); %calc standard dev

if tails
    stats.tails = 'one-sided';
    alpha = conf/100;
    z = tinv(alpha,n-1);
else
    stats.tails = 'two-sided';
    alpha = (conf/100)+((1-(conf/100))/2);
    z = tinv(alpha,n-1);               %choose how many standard dev's for conf interval of 95%
end

stats.conf = conf;
stats.df = n-1;
stats.z = z;

mu = mean(perc_missing); %find mean

if sd_flag
    ci = z*sd;
    stats.method = 'standard deviation';
else
    ci = z*(sd/sqrt(n));    %calculate confidence interval
    stats.method = 'standard error';
end

stats.mu = mu;
stats.upper = mu+ci;
if ~tails
    stats.lower = mu-ci;
end


labels = dso.label{1};
if isempty(labels)
    label_cell = mat2cell(sam_include,ones(size(sam_include,1),1), ones(size(sam_include,2),1));
else
    label_cell = {};
    for i = 1:length(sam_include)
        label_cell{end+1} = labels(sam_include(i),:);
    end
end

h = figure('position', [0, 0, 400, 1600]);
set(h,'Visible','off');

ax = axes;
hold on;

line(ones(1,length(num_missing)+1).*mu,[0,sam_include], 'Color',[0 1 0],'Linestyle','-');
line(ones(1,length(num_missing)+1).*mu+ci,[0,sam_include], 'Color',[0 0 0],'Linestyle','-.');

flag_samples = find(perc_missing > stats.upper);
if ~tails
    line(ones(1,length(num_missing)+1).*mu-ci,[0,sam_include], 'Color',[0 0 0],'Linestyle','-.');
    flag_samples = [flag_samples find(perc_missing < stats.lower)];
end

barh(ax, [1:length(sam_include)], perc_missing);
set(ax, 'YDir', 'reverse');
set(ax, 'YTick', [1:length(sam_include)],'YTickLabel',label_cell, 'TickDir', 'out');
title(ax, 'Percentage of missing values per sample');

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 0.2*length(sam_include)])
print(gcf,'-dpng','-r300',figurefile);

close(h);

if del == 1
    eddata = dso; eddata(flag_samples,:) = []; % Hard delete samples
    autoexport(eddata, outfile, 'xml');
elseif del == 2
    eddata = delsamps(dso,flag_samples);  % Soft delete samples
    autoexport(eddata, outfile, 'xml');
else
    autoexport(dso, outfile, 'xml');
end

end

