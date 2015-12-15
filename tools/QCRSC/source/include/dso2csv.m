function dso2csv(dso_xml,csv_data,csv_peaks)

[dso, ~, ~] = autoimport(dso_xml, 'xml');

[n, m] = size(dso); % Make sure that I ignore blanks here?
data_file.Idx = [1:n]';
count = 2;
for i = 1:n
    if strcmp(dso.classid{1}(i),'QC') || strcmp(dso.classid{1}(i),'qc') || strcmp(dso.classid{1}(i),'Qc') | strcmp(dso.classid{1}(i),'qC')
        data_file.Sample_Type{i} = 'QC';
        data_file.Sample{i} = 0;
        data_file.QC{i} = 1; 
        data_file.Sample_Rep{i} = 1;
    else
        data_file.Sample_Type{i} = 'sample';
        data_file.Sample{i} = 1;
        data_file.QC{i} = 0;
        data_file.Sample_Rep{i} = count;
        count = count + 1;
    end
end
data_file.Sample_Type = data_file.Sample_Type';
data_file.Sample = data_file.Sample';
data_file.QC = data_file.QC';
data_file.Batch = dso.class{1,3}';
data_file.Sample_Rep = data_file.Sample_Rep';
data_file.Order = dso.class{1,4}';
data_file.Class = dso.classid{1}';
data_file.Class_2 = dso.label{1}; %Somehow get a space inbetween class and 2
data_file = orderfields(data_file,{'Idx','Sample_Type','Sample','QC','Batch','Order','Sample_Rep','Class','Class_2'});
field_name = [];
for j = 1:m
    field_name{j} = strcat('M',num2str(j));
    field_name_mz{j} = strcat('M',num2str(j),':',num2str(dso.axisscale{2}(j),'%.5f'));
    data_file.(field_name{j}) = dso.data(:,j); 
end
T = struct2table(data_file);

writetable(T,csv_data,'Delimiter',',');

peak_file.Name = field_name';
peak_file.Label = field_name_mz';
peak_file.Mass = dso.axisscale{2}';
T2 = struct2table(peak_file);
writetable(T2,csv_peaks,'Delimiter',',');
