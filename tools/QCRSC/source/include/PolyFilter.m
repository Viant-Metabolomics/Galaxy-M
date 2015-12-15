function [Data_Filtered,Peaks_Filtered] = PolyFilter(Meta,Data,Peaks,C)


list = Data.Properties.VariableNames;
opts = fitoptions(C.polyFunc, 'Robust', 'Bisquare');


Peaks_Filtered = Peaks;

Xall = table2array(Data);
[rows,cols] = size(Xall);

if (C.zeroflag == 1),  Xall(Xall == 0) = NaN; end;
if strcmp(C.peakTransform,'log'), Xall = log(Xall); end;

Batch_list = unique(Meta.Batch);
Batch_num = numel(Batch_list);


QCcut = nan(cols,Batch_num);
Scut = nan(cols,Batch_num);
Order = Meta.Order;
Mask = zeros(rows,cols);
Mask2 = zeros(rows,cols);
warning ('off','all');
for i = 1:Batch_num
    if C.report, display(['Processing Batch ',num2str(Batch_list(i))]); end;
    batch = find(ismember(Meta.Batch,Batch_list(i)));
    
    X = Xall(batch,:);
    T = Meta.Order(batch);
    QC = Meta.QC(batch);
    
    for j = 1:cols
        %fprintf('%d ',j);
        %if ~rem(j,50), fprintf('\n'); end
        if ~rem(j,100), fprintf('... '); end
        x = X(:,j);
        t = T;
        qc = QC;

        temp = isnan(x);

        x(temp) = [];
        t(temp) = [];
        qc(temp) = [];

        Xsample = x(~qc);
        Tsample = t(~qc);
        
        if length(Tsample) < 5, continue, end; 
        
        mdl = fit(Tsample,Xsample,C.polyFunc,opts);
        CI = predint(mdl,t,C.polyCI);

        Xqc = x(qc);
        Tqc = t(qc);
        CIqc = CI(qc,:);
        CIsample = CI(~qc,:);

        cut_up = Xqc > CIqc(:,2);
        cut_low = Xqc < CIqc(:,1);
        qc_cut = or(cut_up,cut_low);
        Tcut = Tqc(qc_cut);
        cut = ismember(Order,Tcut);
        Mask(cut,j) = 1;
        QCcut(j,i) = sum(qc_cut);

        cut_up = Xsample > CIsample(:,2);
        Tcut = Tsample(cut_up);
        cut = ismember(Order,Tcut);
        Mask2(cut,j) = CIsample(cut_up,2);
        cut_low = Xsample < CIsample(:,1);
        Tcut = Tsample(cut_low);
        cut = ismember(Order,Tcut);
        Mask2(cut,j) = CIsample(cut_low,1);
        Scut(j,i) = sum(or(cut_up,cut_low));


    end
end



warning ('on','all'); 
 

if C.report, display('removing QC outliers'); end;
Mask = logical(Mask);
  
Xall(Mask) = NaN;

M = logical(Mask2);
if strcmp(C.polyAction,'suppress')
    if C.report, display('suppressing sample outliers'); end;
    Xall(M) = Mask2(M);
elseif strcmp(C.polyAction,'remove')
    if C.report, display('removing sample outliers'); end;
    Xall(M) = NaN;
end

if strcmp(C.peakTransform,'log'), Xall = exp(Xall); end;
    
Data_Filtered = array2table(Xall,'VariableNames',list);

for i = 1:Batch_num
    Peaks_Filtered.(['QC_Outliers_B',num2str(Batch_list(i))]) = QCcut(:,i);
    Peaks_Filtered.(['Sample_Outliers_B',num2str(Batch_list(i))]) = Scut(:,i);
end


end