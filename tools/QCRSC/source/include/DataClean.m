function [Data2,Peaks2,Peaks4,T,listx] = DataClean(Meta,Data,Peaks,C)

num_test = 8;

if C.QCRSC == 0
    C.DataFile = ['QCRSC_',C.ProjectName,'_data.csv'];
    C.PeakFile = ['QCRSC_',C.ProjectName,'_peaks.csv'];
    [Meta,Data,Peaks] = load_data(C);
end

Peaks2 = Peaks;

Xall = table2array(Data);

MedMPA = Peaks.MedMPA;

rows = size(Peaks,1);
error_code = zeros(rows,num_test);

Batch_list = unique(Meta.Batch);
Batch_num = numel(Batch_list);
Z = [];
for i = 1:Batch_num
    temp = Peaks.(['MPA_B',num2str(Batch_list(i))]);
    zeta = abs((MedMPA-temp)./MedMPA);
    temp = zeta > C.KsiTol;
    Peaks2.(['Ksi_B',num2str(Batch_list(i))]) = zeta;
    Z = [Z,temp];
    batch = ismember(Meta.Batch,Batch_list(i));
    Xall(batch,temp) = NaN;   
end

list = Data.Properties.VariableNames;
Data2 = array2table(Xall,'VariableNames',list);

Peaks2.MPAkill = sum(Z,2);

temp = Peaks2.MPAkill > 0;

error_code(temp,1) = 1;

temp1 = or(temp,Peaks.SumKILL);

error_code(temp1,1) = 2*10^1;

if Batch_num > 1
    temp2 = Peaks.Krus_p < C.pdiff;
    error_code(temp2,2) = 3*10^2;
else
    temp2 = zeros(numel(temp1),1);
end

if C.kill
    cut = or(temp1,temp2);
else
    cut = temp2;
end

temp3 = Peaks.QC_offset > C.QCdist;

error_code(temp3,3) = 4*10^3;

cut = or(cut,temp3);

%[Peaks2] = RSD_D_RATIO(Meta,Data2,Peaks2,C);

%temp3 = Peaks2.RSD_T > C.MaxRSD;

temp3 = Peaks2.MaxRSD > C.MaxRSD;

error_code(temp3,4) = 5*10^4;

cut = or(cut,temp3);

%temp3 = isnan(Peaks2.RSD_T);
temp3 = isnan(Peaks2.MaxRSD);

error_code(temp3,5) = 6*10^5;

cut = or(cut,temp3);

%temp3 = Peaks2.D_RATIO_T < C.MinD_RATIO;
temp3 = Peaks2.MinD_RATIO < C.MinD_RATIO;

error_code(temp3,6) = 7*10^6;

cut = or(cut,temp3);

%temp3 = isnan(Peaks2.D_RATIO_T);
temp3 = isnan(Peaks2.MinD_RATIO);

error_code(temp3,7) = 8*10^7;

cut = or(cut,temp3);

listx = list(cut);
list = list(~cut);

Xall = Xall(:,~cut);

Data2 = array2table(Xall,'VariableNames',list);

temp = num2str(sum(error_code,2),['%0',num2str(num_test),'i']);
temp11 = char(zeros(size(temp,1),size(temp,2)+2));
for i = 1:size(temp,1)
    temp11(i,:) = ['[',fliplr(temp(i,:)),']']; 
end

Peaks2.cut = double(cut);
Peaks2.ErrorCode = temp11;
Peaks3 = Peaks2(listx,:);
Peaks2 = Peaks2(list,:);

list2 = Peaks3.Properties.VariableNames;
T = table2cell(Peaks3);
 TT = [listx',T];
Peaks4 = cell2table(TT,'variablenames',[{'Name'},list2]);

[pathstr,~,~] = fileparts(C.DataFile);

writetable(Peaks4,[pathstr,filesep,'peaks_removed.csv']);
if C.report, display([num2str(sum(cut)),' Peaks Removed']); end;
end

