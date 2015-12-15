function [Peaks_new] = compare_b2qc(Meta,Data,Peaks,C)

X = table2array(Data);

cols = size(X,2);
if strcmp(C.peakTransform,'log')
    X = log(X);
end

QC = Meta.QC;

QC_offset = nan(cols,1);

fprintf('.');
for i = 1:cols
    if ~rem(i,100), fprintf('.'); end
    
    medQC = nanmedian(X(QC,i));
    medSample = nanmedian(X(~QC,i));
    diff = abs(medSample-medQC);
    madSample = mad(X(~QC,i));
    QC_offset(i) = diff/madSample;
end
fprintf('\n');
Peaks_new = Peaks;
Peaks_new.QC_offset = QC_offset;

end

