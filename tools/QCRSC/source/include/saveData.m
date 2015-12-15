function [] = saveData(Meta,Data,Peak,name)

T = table2cell(Data);
temp = cell2mat(T);
temp(isnan(temp)) = -99;
T = num2cell(temp);
T2 = table2cell(Meta);
TT = [T2,T];
list = Data.Properties.VariableNames;
list2 = Meta.Properties.VariableNames;
Meta_Data = cell2table(TT,'variablenames',[list2,list]);
writetable(Meta_Data,[name,'_data.csv']);

list = Peak.Properties.RowNames;
list2 = Peak.Properties.VariableNames;
T = table2cell(Peak);
TT = [list,T];
Peaks = cell2table(TT,'variablenames',[{'Name'},list2]);
writetable(Peaks,[name,'_peaks.csv']);

       
end

