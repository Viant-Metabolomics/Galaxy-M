function [Peaks] = RSD_SNR(Meta,Data,Peaks,C)


X = table2array(Data);
QC = Meta.QC;

X1qc = X(QC,:);
X1sample = X(~QC,:);

if strcmp(C.peakTransform,'log')

    MPAqc = nanmedian(X1qc);
    RSD = 1.4826*100*mad(X1qc,1)./MPAqc;
    SNR = 20*log10(mad(X1sample,1)./mad(X1qc,1));
    
    temp3 = ~isnan(X1qc);
    numQCs = sum(temp3);

    % This is a cludge for when the number of QCs is so low that MAD becomes
    % meaningless


    if numQCs < 5
        RSD = 100*std(X1qc)/nanmean(X1qc);
        SNR = 20*log10(nanstd(X1sample)/nanstd(X1qc));
    end
    
else
    MPAqc = nanmedian(X1qc);
    RSD = 100*std(X1qc)./nanmean(X1qc);
    SNR = 20*log10(mad(X1sample)./mad(X1qc));
end

Peaks.MPA_T = MPAqc';
Peaks.RSD_T = RSD';
Peaks.SNR_T = SNR';

end

