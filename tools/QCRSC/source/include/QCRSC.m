function [Z,Report,xx] = QCRSC(X,T,QC,C)

params = C.searchRange(1):C.searchRange(2):C.searchRange(3);
QCwindow = C.QCwindow;
max_window = C.maxWindow;

kill = 0;

TTqc = T(QC);
XXqc = X(QC);

Xqc = XXqc(~isnan(XXqc));
Tqc = TTqc(~isnan(XXqc));
MPA = median(Xqc);
numQCs = numel(Tqc);

if numQCs > 2   
    
    leadin = floor((Tqc(1) - TTqc(1))/QCwindow);
    leadout = floor((TTqc(end) - Tqc(end))/QCwindow);

    if leadin > 2
        Xqc = [Xqc(1);Xqc];
        Tqc = [TTqc(1);Tqc];
    end

    if leadout > 2
        Xqc = [Xqc;Xqc(end)];
        Tqc = [Tqc;TTqc(end)];
    end
    
    temp2 = [Tqc(2:end);Tqc(end)];
    maxdist = nanmax(temp2-Tqc)-1;
    MAXwindow = floor((maxdist/(QCwindow))-1);
else
    MAXwindow = numel(TTqc);
    leadin = numel(TTqc);
    leadout = numel(TTqc);
    maxdist = numel(TTqc);
end


if leadin < max_window && leadout < max_window
    h = ((Tqc(end)-Tqc(1))/(numel(Tqc)-1));
    epsilon = h^3/16;
    if length(Xqc) < 5 
        % QCs < 5 cannot effectively perform QCspline cross-valiadation
        % setting opt_param to effectively a linear correction.
        opt_param = max(params);
        cvMse = nan(1,numel(params));
        Report.Poly = 0;
    else
        Report.Poly = 1;
        cvMse = nan(numel(params),1);
        for i = 1:numel(params)
            p = 1/(1+epsilon*10^(params(i)));
            regfr=@(Ttrain,Xtrain,Ttest)(csaps(Ttrain,Xtrain,p,Ttest));
            cp = cvpartition(length(Xqc),'leaveout');
            cvMse(i) = crossval('mse',Tqc,Xqc,'predfun',regfr,'partition',cp);
        end
        [dummy,min_cvMse] = min(cvMse);
        opt_param = params(min_cvMse);
    end
    
    p = 1/(1+epsilon*10^(opt_param));
    Z = csaps(Tqc,Xqc,p,T);
    
else
    % There are too few QCs so peak being replaced with missing values
    Report.Poly = -1;
    Z = nan(length(X),1);
    opt_param = NaN;
    MPA = NaN;
    cvMse = nan(1,numel(params));
    kill = 1;
end    

if strcmp(C.operator,'divide')
    zz = X./Z;
    xx = zz.*MPA;
else
   zz = X-Z;
   xx = zz+MPA;
end

if strcmp(C.peakTransform,'log')
    xx = exp(xx);
end;

    xxqc = xx(QC);
    xxsample = xx(~QC);
    
    
    
%      RSD = 100*nanstd(xxqc)/nanmean(xxqc);
%      SNR = 20*log10(nanstd(xxsample)/nanstd(xxqc));
      RSD = 100*1.4826*mad(xxqc,1)/nanmedian(xxqc);
      D_RATIO = mad(xxsample,1)/mad(xxqc,1);
      
      
% This is a cludge for when the number of QCs is so low that MAD becomes
% meaningless


if numQCs < 5
    RSD = 100*nanstd(xxqc)/nanmean(xxqc);
    D_RATIO = nanstd(xxsample)/nanstd(xxqc);
end
      
      
      if ~isfinite(D_RATIO)
           D_RATIO = NaN;
      end

Report.kill = kill;
Report.QCwindow = QCwindow;
Report.cvMse = cvMse; 
Report.opt_param = opt_param;
Report.MPA = MPA;
Report.numQCs = numQCs;
Report.leadin = leadin;
Report.leadout = leadout;
Report.maxdist = maxdist;
Report.MAXwindow = MAXwindow;
Report.RSD = RSD;
Report.D_RATIO = D_RATIO;
end

