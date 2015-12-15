function [Peaks_new] = compare_batches(Meta,Data,Peaks,C)

X = table2array(Data);

if strcmp(C.peakTransform,'log')
    X = log(X);
end

isQC = Meta.QC;

cols = size(X,2);
Krus_p = nan(cols,1);
ANOVA_p = nan(cols,1);

display('Comparing Batches Statistically');
fprintf('.');
for i = 1:cols
    if ~rem(i,100), fprintf('.'); end
    if Peaks.SumKILL(i) > 0, continue; end;
    Krus_p(i) = kruskalwallis(X(~isQC,i),Meta.Batch(~isQC),'off');
    ANOVA_p(i) = anova1(X(~isQC,i),Meta.Batch(~isQC),'off');
end
fprintf('\n');

Peaks_new = Peaks;
Peaks_new.Krus_p = Krus_p;
Peaks_new.ANOVA_p = ANOVA_p;

end

