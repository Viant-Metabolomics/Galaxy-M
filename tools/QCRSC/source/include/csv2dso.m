function csv2dso(csv_data,csv_peaks,dso_xml)

%% Load peak file
Peaks = readtable(csv_peaks,'ReadRowNames',true);

%% Load data file
DataTable = readtable(csv_data,'TreatAsEmpty',{'','.','NA','N/A'});
data_matrix = table2array(DataTable(:,10:end));

%% Create data set object
dso = dataset(data_matrix);

dso.axisscale{2} = table2array(Peaks(:,2));

dso.label{1} =  table2array(DataTable(:,9));
dso.labelname{1} = 'sample_name';

dso.classname{1,1} = 'classes';
dso.classid{1,1} = table2array(DataTable(:,8));

dso.classname{1,2} = 'blank_flag';
classnames_BF ={};
for i = 1:size(data_matrix,1);
    if strncmpi(dso.classid{1,1}(i),'blank',length('blank'));
        classnames_BF{end+1} = 'blank';
    else
        classnames_BF{end+1} = 'sample';
    end
end
dso.classid{1,2} = classnames_BF;

dso.classname{1,3} = 'batch_info';
dso.classid{1,3} = table2array(DataTable(:,5));

dso.classname{1,4} = 'order_info'; 
dso.classid{1,4} = table2array(DataTable(:,6));

autoexport(dso, dso_xml, 'xml');
