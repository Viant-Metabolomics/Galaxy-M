function [Data_Filtered,Peaks_Filtered] = PolyFilter_ALL(Meta,Data,Peaks,C)

Xall = table2array(Data);
[rows,cols] = size(Xall);

if C.zeroflag == 1
    Xall(Xall == 0) = NaN;
end

if strcmp(C.peakTransform,'log')
    Xall(Xall == 0) = NaN;
    Xall = log(Xall); 
end;


Batch_list = unique(Meta.Batch);
Batch_num = numel(Batch_list);
warning ('off','all');


for i = 1:Batch_num
    if C.report, display(['Processing Batch ',num2str(Batch_list(i))]); end;
    batch = find(ismember(Data.Batch,Batch_list(i)));
    
    X = Xall(batch,:);
    T = Meta.Order(batch);
    QC = Meta.QC(batch);
    if parallel_flag
        if C.report, display('reporting is disabled due to parallel processing ... be very patient'); end 
        parfor j = 1:cols
            [Xout(:,j),R(j)] = PolyFilter(X(:,j),T,QC,C);
        end
    else
        if C.report, display('Maybe you should buy the Parallel Processing Toolbox to speed things up!'); end;
        for j = 1:cols
            if C.report, if ~rem(j,100), fprintf('...'); end; end;
            [Xout(:,j),R(j)] = PolyFilter(X(:,j),T,QC,C);
        end
    end
    if C.report, fprintf('\n'); end;
    
end
                
list = Data.Properties.VariableNames;              
Data_Filtered = array2table(Xout,'VariableNames',list);

Poly = struct2table(R,'RowNames',list);
Peaks_Filtered = join(Peaks,Poly,'Keys','RowNames');

end