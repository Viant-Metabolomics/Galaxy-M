function [C] = load_config(config_file)
%Parse the QCRSC configuration file into a Struct
%With validation
%   config_file = configuration text file
%   here is an example:
% -------------------------------------------------------------------------
% ProjectName="Name for the project - used in all output"(txt)=Test1
% DataFile=".csv filename containing the raw data"(*.csv)=Data.csv
% PeakFile=".csv filename for the associated peak data"(*.csv)=Peak.csv
% PolyFilter="Do you want to perform Prefilter?"(yes/no)=yes
% QCRSC="Do you want to perform QCRSC?"(yes/no)=yes
% DataClean="Do you want to perform Data cleaning?"(yes/no)=yes
% 
% %%%%%%%%%%% General Options %%%%%%%%%%%%%%%%%%%%
% 
% report="verbose reporting"(yes/no)=yes
% parallel="using Matlab parallel processing?"(yes/no)=yes
% saveopt="save at each stage"(yes/no)=yes
% zeroflag="replace zeros with NaNs"(yes/no)=yes
% peakTransform="peak transform before any processing"(log/none)=log
% plotflag="automatically save figures as tiff files"(yes/no)=yes
% 
% %%%%%%%%%%%% PolyFilter options %%%%%%%%%%%%%%%%
% 
% polyFunc="order of polynomial fit"(1:6)=poly3
% polyCI="percentage confidence interval"(double)=0.99
% polyAction="what to do with NON-QC outliers"(suppress/remove/ignore) = suppress
% 
% %%%%%%%%%%%% QCRSC options %%%%%%%%%%%%%%%%
% 
% operator="correction type"(subtract/divide)=subtract
% maxWindow="maximum number of consecutive missing QCs"(int)=-1
% searchRange="search range for B-spline smoothing parameter"(min:step:max)=0:0.5:7
% 
% 
% %%%%%%%%%%%% Data Clean options %%%%%%%%%%%%%%%%
% 
% MPAtol="tolerance for between MPA comparison"(num)=4
% kill="remove peaks with missing batches"(yes/no)=yes
% pdiff="critical p-value for batch comparison"(number)=1e-4
% QCdist="maximum allowable QC offset - number of SDs"(number)=3
% MaxRSD="maximum allowed %RSD"(number)=25
% MinD_RATIO="minimum allowed Dispersion Ratio"(number)=1
%
%
% -------------------------------------------------------------------------
% Can remove the discription if desired e.g.
%
%     peakTransform==log
%     operator=subtract
%     maxWindow=10
%     searchRange=-1:0.5:6
%     report=yes
%


% default settings

C.ProjectName = '';
C.DataFile = '';
C.PeakFile = '';
C.PolyFilter = 0;
C.QCRSC = 1;
C.DataClean = 0;

C.report = 1;
C.parallel = 0;
C.saveopt = 1;
C.zeroflag = 1;
C.peakTransform = 'log';
C.plotflag = 1;

C.polyFunc = 'poly3';
C.polyCI = 99;
C.polyAction = 'ignore';

C.operator = 'subtract';
C.maxWindow = -1;
C.searchRange = [0,0.5,7];


C.MPAtol = 4;
C.kill = 1;
C.pdiff = 1e-4;
C.QCdist = 3;
C.MaxRSD = 25;
C.MinD_RATIO = 0;


if ~isdeployed
    % Check version & toolbox
    v = ver;
    toolboxes = {v.Name};
    if ~any(ismember(toolboxes,'Statistics Toolbox'))
        error('You do not have the Statistics Toolbox');
    end

    if ~any(ismember(toolboxes,'Curve Fitting Toolbox'))
        error('You do not have the Curve Fitting Toolbox');
    end
end


% Parse Config file
if ~ischar(config_file), error('input has to be a string - i.e. ''config.txt'' '); end;
if ~exist(config_file,'file'), error([config_file,' does not exist']); end;
fid = fopen(config_file);

i = 1;

try
    while ~feof(fid)
        ss = fgetl(fid);
        if isempty(ss), continue; end;
        if strcmp(ss(1),'%'), continue; end;
        T(i,:) = textscan(ss, '%s %s %s','Delimiter','=');
        i = i+1;
    end
catch
    error(['configuration file is not in the correct format: line ',num2str(i)]);
end

cols = size(T,2);

if cols == 3
    data = [[T{:,1}]',[T{:,3}]'];
elseif cols == 2
    data = [[T{:,1}]',[T{:,2}]'];
else
    error('configuration file is not in the correct format');
end


for i = 1:size(data,1)
    switch data{i,1}
        case 'ProjectName'
            val = strtrim(data{i,2});
            C.ProjectName = val;
        case 'DataFile'
            val = strtrim(data{i,2});
            C.DataFile = val;
            if ~exist(val,'file'), error([val,' does not exist']); end;
        case 'PeakFile'
            val = strtrim(data{i,2});
            C.PeakFile = val;
            if ~exist(val,'file'), error([val,' does not exist']); end;
        case 'PolyFilter'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.PolyFilter = 1;
            elseif strcmpi(val,'no');
                C.PolyFilter = 0;
            else
                error('illegal value: Prefilter = yes or no');
            end
        case 'QCRSC'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.QCRSC = 1;
            elseif strcmpi(val,'no');
                C.QCRSC = 0;
            else
                error('illegal value: QCRSC = yes or no');
            end
        case 'DataClean'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.DataClean = 1;
            elseif strcmpi(val,'no');
                C.DataClean = 0;
            else
                error('illegal value: DataClean = yes or no');
            end
        case 'report'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.report = 1;
            elseif strcmpi(val,'no');
                C.report = 0;
            else
                error('illegal value: report = yes or no');
            end
        case 'parallel'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.parallel = 1;
                if ~any(ismember(toolboxes,'Parallel Computing Toolbox'))
                    display('You do not have the Parallel Computing Toolbox');
                    C.parallel = 0;
                end
            elseif strcmpi(val,'no');
                C.parallel = 0;
            else
                error('illegal value: parallel = yes or no');
            end     
        case 'saveopt'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.saveopt = 1;
            elseif strcmpi(val,'no');
                C.saveopt = 0;
            else
                error('illegal value: saveopt = yes or no');
            end
        case 'zeroflag'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.zeroflag = 1;
            elseif strcmpi(val,'no');
                C.zeroflag = 0;
            else
                error('illegal value: zeroflag = yes or no');
            end
        case 'peakTransform'
            val = strtrim(data{i,2});
            legal = {'log','none'};
            if sum(ismember(legal,val)) == 0
                error('illegal value: peak_transform = log/none');
            end
            C.peakTransform = val;
        case 'plotflag'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.plotflag = 1;
            elseif strcmpi(val,'no');
                C.plotflag = 0;
            else
                error('illegal value: plotflag = yes or no');
            end
        
        case 'polyFunc'
            val = strtrim(data{i,2});
            legal = {'poly1','poly2','poly3','poly4','poly5','poly6'}; 
            if sum(ismember(legal,val)) == 0
                error('illegal value: poly_func = poly1 , poly2 ... poly6');
            end
            C.polyFunc = val;
        case 'polyCI'
            val = strtrim(data{i,2});
            val = str2double(val);
            if isempty(val)
                error('illegal value: ci must be a number');
            end
            if val > 1
                error('illegal value: ci must be less than 1');
            end
            if val < 0
                error('illegal value: ci must be greater than 0');
            end
            if val < 0.9
                warning('ci is less than 0.9. ARE YOU SURE?');
            end
            C.polyCI = val;   
        case 'polyAction'
            val = strtrim(data{i,2});
            legal = {'suppress','remove','ignore'}; 
            if sum(ismember(legal,val)) == 0
                error('illegal value: action = suppress/remove/ignore');
            end
            C.ployAction = val;   
        case 'operator'
            val = strtrim(data{i,2});
            legal = {'divide','subtract'};
            if sum(ismember(legal,val)) == 0
                error('illegal value: operator = divide/subtract');
            end
            C.operator = val;
        case 'maxWindow'
            val = strtrim(data{i,2});
            val = round(str2double(val));
            if isnan(val)
                error('illegal value: maxdist must be an integer');
            end
            if val <  -1
                error('illegal value: maxdist must be a positive integer');
            end
            if val > 20
                warning('maxdist is greater than 20 ARE YOU SURE?');
            end
            C.maxWindow = val;
        case 'searchRange'
            val = sscanf(strtrim(data{i,2}),'%f:%f:%f');
            if length(val) < 3 || val(3) < val(1) || val(2) < 0 
                error('illegal value: search_range must be in the format start:interval:end e.g. 1:0.5:10');
            end
            if val(2) < 0.1
                warning('search_range step < 0.1');
            end
            if val(2) > val(3)
                warning('search_range step seems very large. ARE YOU SURE?');
            end
            C.searchRange = val;
        case 'KsiTol'
            val = strtrim(data{i,2});
            val = str2double(val);
            if isnan(val)
                error('illegal value: KsiTol must be an real positive number');
            end
            if val <  1
                error('illegal value: KsiTol must be a number > 1');
            end
            if val > 10
                warning('KsiTol is greater than 20 ARE YOU SURE?');
            end
            C.KsiTol = val;
        case 'kill'
            val = strtrim(data{i,2});
            if strcmpi(val,'yes');
                C.kill = 1;
            elseif strcmpi(val,'no');
                C.kill = 0;
            else
                error('illegal value: kill = yes or no');
            end
        case 'pdiff'
            val = strtrim(data{i,2});
            val = str2double(val);
            if isnan(val)
                error('illegal value: pdiff must be a real positive number');
            end
            if val >  1
                error('illegal value: pdiff must be < 1 - it is a p-value!!');
            end
            if val > 0.05
                warning('pdiff is greater than 0.05 ARE YOU SURE?');
            end
            C.pdiff = val;
        case 'QCdist'
            val = strtrim(data{i,2});
            val = str2double(val);
            if isnan(val)
                error('illegal value: QCdist must be an real positive number');
            end
            if val < 0
                error('illegal value: QCdist must be a real positive number');
            end
            if val > 5
                warning('QCdist is greater than 5 ARE YOU SURE?');
            end
            C.QCdist = val;
        case 'MaxRSD'
            val = strtrim(data{i,2});
            val = str2double(val);
            if isnan(val)
                error('illegal value: MaxRSD must be an real positive number');
            end
            if val < 0
                error('illegal value: MaxRSD must be a real positive number');
            end
            if val > 50
                warning('MaxRSD is greater than 50 ARE YOU SURE?');
            end
            C.MaxRSD = val;
        case 'MinD_RATIO'
            val = strtrim(data{i,2});
            val = str2double(val);
            if isnan(val)
                error('illegal value: MinSNR must be an real positive number');
            end
            if val < 0
                error('illegal value: MinSNR must be a real positive number');
            end
            if val > 5
                warning('MinSNR is greater than 5 ARE YOU SURE?');
            end
            C.MinD_RATIO = val;    
        case ''
            continue;
        otherwise
            error(['illegal configuration option: ' data{i,1}]);
    end

end           

if isempty(C.ProjectName)
    error('You must supply a ProjectName');
end
if isempty(C.DataFile)
    error('You must supply a DataFile name');
end
if isempty(C.PeakFile)
    error('You must supply a PeakFile name');
end

temp = C.PolyFilter + C.QCRSC + C.DataClean;

if temp == 0, error('You must perform at least one function: Prefilter, QCRSC and/or DataClean'); end;

ff = fieldnames(C);
if C.report
    
    display('Loading Configuration file ...');
    display(' ');
    
    for i = 1:length(ff)
        if isa(C.(ff{i}),'numeric')
            if strcmp(ff{i},'searchRange')
                display([ff{i},' = ',num2str(C.(ff{i})(1)),':',num2str(C.(ff{i})(2)),':',num2str(C.(ff{i})(3))]);
            else   
                display([ff{i},' = ',num2str(C.(ff{i}))]);
            end
        else
            display([ff{i},' = ',C.(ff{i})]);
        end
    end
    
    display(' ');
    display('... done.');
    display(' ');
end
    
end

