function [Data_QCRSC,Peaks_QCRSC] = QCRSC_ALL(Meta,Data,Peaks,C)

list = Data.Properties.VariableNames;
name = C.ProjectName;
parallel_flag = C.parallel;

Xall = table2array(Data);
[rows,cols] = size(Xall);

if C.zeroflag == 1
    Xall(Xall == 0) = NaN;
end

if strcmp(C.peakTransform,'log')
    Xall(Xall == 0) = NaN;
    Xall = log(Xall); 
end;

Z = nan(size(Xall));


Batch_list = unique(Meta.Batch);
Batch_num = numel(Batch_list);



% MPA = nan(cols,Batch_num);
% OPTPARAM = nan(cols,Batch_num);
% MAXWINDOW = nan(cols,Batch_num);
% RSD = nan(cols,Batch_num);
% D_RATIO = nan(cols,Batch_num);
% KILL = nan(cols,Batch_num);

for i = 1:Batch_num
    if C.report, display(['Processing Batch ',num2str(Batch_list(i))]); end;
    batch = find(ismember(Meta.Batch,Batch_list(i)));    
    X = Xall(batch,:);
    T = Meta.Order(batch);
    QC = Meta.QC(batch);
    
    Tqc = T(QC);
    temp = [Tqc(2:end);Tqc(end)];
    C.QCwindow = median((temp-Tqc)-1);
    
    if C.maxWindow == -1,
        C.maxWindow = floor(numel(Tqc)/2);
    end
    
    if C.report, display(['max_window ',num2str(C.maxWindow)]); end;
    if C.report, display('Correcting peaks ... this make take a long time\n'); end 
    
    if parallel_flag
        if C.report, display('reporting is disabled due to parallel processing ... be very patient'); end 
        parfor j = 1:cols
            [Z(batch,j),RR(j)] = QCRSC(X(:,j),T,QC,C); 
        end
    else
        if C.report, display('Maybe you should buy the Parallel Processing Toolbox to speed things up!'); end;
        for j = 1:cols
            if C.report, if ~rem(j,100), fprintf('...'); end; end;
            [Z(batch,j),RR(j)] = QCRSC(X(:,j),T,QC,C); 
        end
    end
    if C.report, fprintf('\n'); end;
    RSC{i} = struct2table(RR,'RowNames',list);
    Batch_QCRSC = join(Peaks,RSC{i},'Keys','RowNames');
    
    [pathstr,~,~] = fileparts(C.DataFile);
    filename = [pathstr,filesep,name,'_Batch_',num2str(Batch_list(i)),'_Report.csv'];
    writetable(Batch_QCRSC,filename);
end



for i = 1:Batch_num
        TMPA(:,i) = RSC{i}.MPA;
        TRSD(:,i) = RSC{i}.RSD;
        TD_RATIO(:,i) = RSC{i}.D_RATIO;
        TKILL(:,i) = RSC{i}.kill;
        TMAXWIN(:,i) = RSC{i}.MAXwindow;
        TOPTPARAM(:,i) = RSC{i}.opt_param;
end


for i = 1:cols
KILLSTR{i,1} = ['%[',num2str(TKILL(i,:)),']'];
end
R.KILLSTR = KILLSTR;
R.SumKILL = sum(TKILL,2);
R.MaxWin = nanmax(TMAXWIN,[],2);
R.MedParam = nanmedian(TOPTPARAM,2);
R.MedMPA = nanmedian(TMPA,2);
R.MaxRSD = nanmax(TRSD,[],2);
R.MinD_RATIO = nanmin(TD_RATIO,[],2);

if strcmp(C.operator,'divide')
    XX = Xall./Z;
    XX = XX.*(ones(rows,1)*R.MedMPA');
else
   XX = Xall-Z;
   XX = XX+(ones(rows,1)*R.MedMPA');
end

if strcmp(C.peakTransform,'log'), 
    XX = exp(XX);
    R.MedMPA = exp(R.MedMPA);
    TMPA = exp(TMPA);
end;




Data_QCRSC = array2table(XX,'VariableNames',list);

for i = 1:Batch_num
    R.(['MPA_B',num2str(Batch_list(i))]) = TMPA(:,i);
end
for i = 1:Batch_num
    R.(['OPTPARAM_B',num2str(Batch_list(i))]) = TOPTPARAM(:,i);
end
for i = 1:Batch_num
    R.(['MAXwin_B',num2str(Batch_list(i))]) = TMAXWIN(:,i);
end
for i = 1:Batch_num
    R.(['RSD_B',num2str(Batch_list(i))]) = TRSD(:,i);
end
for i = 1:Batch_num
    R.(['D_RATIO_B',num2str(Batch_list(i))]) = TD_RATIO(:,i);
end

if Batch_num > 1
    [R] = compare_batches(Meta,Data_QCRSC,R,C);
end

[R] = compare_b2qc(Meta,Data_QCRSC,R,C);

RR = struct2table(R,'RowNames',list);
Peaks_QCRSC = join(Peaks,RR,'Keys','RowNames');

% plot_kds(D_RATIO,Batch_list,'D_RATIO',1);
% plot_kds(RSD,Batch_list,'RSD',1);
% plot_kds(MPA,Batch_list,'MPA',1);
    


end
