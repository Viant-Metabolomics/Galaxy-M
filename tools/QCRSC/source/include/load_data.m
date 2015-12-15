function [Meta,Data,Peaks] = load_data(C)

% This function loads and validates the DataFile and PeakFile DataFile:
% IMPORTANTLY it also sorts the data by ORDER. 
% The first 4 columns must contain metadata - Idx ; Order ; QC ; Batch 
% Metabolite IDs must start with 'M' ... best to use M1 M2 M3 M4 etc. 
% Remaining columns are assumed to be user specific meta data and are ignored. 
% Peak File: The first columns should contain the Peak Label matching the DataFile (M1 M2 .. )
% The remaining columns can contain anything you like. Statistics will be added to this "table"
%
% zeroflag = 1 -> replace zeros with Nans


% Import Data

filename1 = C.DataFile;
filename2 = C.PeakFile;

if exist(filename1,'file') ~= 2
    error(['Filename: ',filename1,' does not exist']);
end
if exist(filename2,'file') ~= 2
    error(['Filename: ',filename2,' does not exist']);
end

[pathstr,name,ext] = fileparts(filename1); 
if ~strcmp(ext,'.csv')
    error(['Filename: ',filename1,' should be a .csv file']);
end
[pathstr,name,ext] = fileparts(filename2); 
if ~strcmp(ext,'.csv')
    error(['Filename: ',filename2,' should be a .csv file']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD PEAK FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if C.report
    display(['Loading PeakFile: ',filename2]);
end

Peaks = readtable(filename2,'ReadRowNames',true);
listP = Peaks.Properties.RowNames;


peaks = unique(listP);
if numel(peaks) ~= numel(listP)
    error(['All Peak Names in ',filename2',' should be unique and match the Peak Names in ',filename2]);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if C.report
    display(['Loading DataFile: ',filename1]);
end

DataTable = readtable(filename1,'TreatAsEmpty',{'','.','NA','N/A'});
list = DataTable.Properties.VariableNames;
found = find(ismember(list,'Order'));
if isempty(found), error('Data file does not contain the required "Order" column'); end;
if sum(isnan(DataTable.Order)), error('Order column cannot contain missing values'); end;
if numel(unique(DataTable.Order)) ~= numel(DataTable.Order)
    error('Order values are not unique. Please change');
end

% IMPORTANT ::::: Sorting all data by INJECTION ORDER

DataTable = sortrows(DataTable,'Order');

found = find(ismember(list,'Idx'));
if isempty(found), error('Data file does not contain the required "Idx" column'); end;
if sum(isnan(DataTable.Idx)), error('Idx column cannot contain missing values'); end;
if numel(unique(DataTable.Idx)) ~= numel(DataTable.Idx)
    error('Idx numbers are not unique. Please change');
end

found = find(ismember(list,'Batch'));
if isempty(found), error('Data file does not contain the required "Batch" column'); end;
if sum(isnan(DataTable.Batch)), error('Batch column cannot contain missing values'); end;

found = find(ismember(list,'QC'));
if isempty(found), error('Data file does not contain the required "QC" column'); end;
if sum(isnan(DataTable.QC)), error('QC column cannot contain missing values'); end;

temp = unique(DataTable.QC);

if numel(temp)~=2
    error('QC values must be binary (0 or 1) & no missing values');
end

if temp(1) ~= 0 && temp(2) ~= 1, error('QC values must be binary (0 or 1) & no missing values'); end;

DataTable.QC = logical(DataTable.QC);

temp = ismember(listP,list);

if sum(temp) ~= numel(listP)
    error(['The Peak Names in ',filename1,' should exactly match the Peak Names in ',filename2,' ( M1, M2 etc. ) Remember that all Peak Names should be unique']);
end


temp = ismember(list,listP);
listA = list(~temp);
listB = list(temp);

try
    X = DataTable{:,listB};
catch
     error(['Fatal Error: Could not load complete data_file.' ... 
        '. You probably have text in the DataFile (other than the header) - ', ...
        'this is not allowed all data must be numberic.', ...
        'Remember the legal form of missing values is either ,, or ,NaN, Please check your file.']);
end

Data = DataTable(:,listB);
Data = standardizeMissing(Data,{NaN,-99,'','.','NA','N/A'});
Meta = DataTable(:,listA);

if C.report
    b = numel(unique(Meta.Batch));
    [r,c] = size(X);
    display(['TOTAL SAMPLES: ',num2str(r),' TOTAL PEAKS: ',num2str(c),' NUMBER OF BATCHES: ',num2str(b)]);
end




if C.report
    display('Done!');
end

end

