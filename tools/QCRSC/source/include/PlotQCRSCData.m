function [fns_jpeg] = PlotQCRSCData(config_filename, peaks_to_plot, outdir)

[C] = load_config(config_filename);
[Meta,Data,Peaks] = load_data(C);
list = Peaks.Properties.RowNames;

Batchlist = unique(Meta.Batch);
fns_jpeg = cell(1,length(peaks_to_plot) * length(Batchlist));
fign = 0;
for k = 1:length(peaks_to_plot)

    if ~isa(peaks_to_plot(k), 'char')
        peak = num2str(peaks_to_plot(k));
    else
        peak = peaks_to_plot(k);
    end
    
    if isempty(strfind(peak, 'M'))
        peak = strcat('M',peak);
    end;
    index = ismember(list,peak);
    if sum(index) == 0
        display(['Fatal Error: peak ',peak,' does not exist']);
        continue;
    end;
    if sum(index) > 1
        display(['Fatal Error: more than one peak called',peak]);
        continue;
    end;
    for j = 1:length(Batchlist)
        
        fign = fign + 1;
        
        batch_num = Batchlist(j);
        batch = ismember(Meta.Batch,batch_num);
        
        X = Data.(peak);
        X = X(batch);
        X2 = X;
        T = Meta.Order(batch);
        QC = Meta.QC(batch);
        if (C.zeroflag == 1),  X(X == 0) = NaN; end;
        if strcmp(C.peakTransform,'log'),  X = log(X); end;
        
        % calculate QC window
        Tqc = T(QC);
        temp = [Tqc(2:end);Tqc(end)];
        C.QCwindow = mean((temp-Tqc)-1);
        
        if C.maxWindow == -1, C.maxWindow = floor(numel(Tqc)/2); end
        
        params = C.searchRange(1):C.searchRange(2):C.searchRange(3);
        
        [z,R,XX] = QCRSC(X,T,QC,C);
        
        Xqc = X2(QC);
        Xsample = X2(~QC);
        
        XXqc = XX(QC);
        XXsample = XX(~QC);
        nanstd(XXsample);
        nanstd(XXqc);
                
        R.RSD = 100*mad(XXqc,1)/nanmedian(XXqc);
        R.SNR = 20*log10(mad(XXsample,1)/mad(XXqc,1));
        BRSD = 100*1.4826*mad(Xqc,1)/nanmedian(Xqc);
        BSNR = 20*log10(mad(Xsample,1)/mad(Xqc,1));
        
        if R.numQCs < 5
            R.RSD = 100*nanstd(XXqc)/nanmean(XXqc);
            R.SNR = 20*log10(nanstd(XXsample)/nanstd(XXqc));
            BRSD = 100*nanstd(Xqc)/nanmean(Xqc);
            BSNR = 20*log10(nanstd(Xsample)/nanstd(Xqc));
        end
        
        %---------------------------------------------------------------
        col = [255 10 71; 51 204 255]./255;
        Label1 = cell(size(QC));
        for i = 1:numel(QC)
            if QC(i) == 1, Label1{i} = 'QC'; else Label1{i} = 'Sample'; end
        end
        
        scrsz = get(0,'ScreenSize');
        
        h = figure('Position',[(scrsz(3)*1/40) scrsz(4)/8 scrsz(3)*0.6 scrsz(4)*0.4]);
        set(h,'Visible','off');
        subplot(2,3,2:3)
        if strcmp(C.peakTransform,'log')
            hold on;
            hqc = scatter(T(QC),X(QC),100,'fill','o','MarkerFaceColor',col(2,:),'MarkerEdgeColor','k');
            hs = scatter(T(~QC),X(~QC),100,'fill','s','MarkerFaceColor',col(1,:),'MarkerEdgeColor','k');
            hold off;
            ylabel('log(Peak Area)');
        else
            hold on;
            hqc = scatter(T(QC),X(QC),100,'fill','o','MarkerFaceColor',col(2,:),'MarkerEdgeColor','k');
            hs = scatter(T(~QC),X(~QC),100,'fill','s','MarkerFaceColor',col(1,:),'MarkerEdgeColor','k');
            hold off;
            ylabel('Peak Area');
        end
        
        title(['Batch ',num2str(batch_num),'; Peak ',peak,'; m/z ',num2str(table2array(Peaks(peak,'Mass')),'%.5f'),' -- BEFORE CORRECTION']);
        hold on;
        
        r = max(T);
        qcmean = nanmean(X(QC));
        qcstd = nanstd(X(QC));
        samplemean = nanmean(X(~QC));
        samplestd = nanstd(X(~QC));
        
        axis([min(T)-1, max(T)+1, min([nanmin(X),samplemean-(2*samplestd)])*0.99, max([nanmax(X),samplemean+(2*samplestd)])*1.01]);
        
        plot([0,r+1]',[samplemean+(2*samplestd),samplemean+(2*samplestd)],'--k');
        text(r+1,samplemean+(2*samplestd),'+2SD (Sample)');
        plot([0,r+1]',[samplemean-(2*samplestd),samplemean-(2*samplestd)],'--k');
        text(r+1,samplemean-(2*samplestd),'-2SD (Sample)');
        plot([0,r+1]',[qcmean+(2*qcstd),qcmean+(2*qcstd)],'--b');
        text(r+1,qcmean+(2*qcstd),'+2SD (QC)');
        plot([0,r+1]',[qcmean-(2*qcstd),qcmean-(2*qcstd)],'--b');
        text(r+1,qcmean-(2*qcstd),'-2SD (QC)');
        plot([0,r+1]',[qcmean,qcmean],'-.b');
        text(r+1,qcmean,'Mean QC');
        plot([0,r+1]',[samplemean,samplemean],'-.k');
        text(r+1,samplemean,'Mean Sample');
                
        plot(T,z,'k^')
        hold off;

        if nansum(XX) ~= 0
                      
            subplot(2,3,5:6)
            
            if strcmp(C.peakTransform,'log')
                XX = log(XX);
                hold on;
                hqc = scatter(T(QC),XX(QC),100,'fill','o','MarkerFaceColor',col(2,:),'MarkerEdgeColor','k');
                hs = scatter(T(~QC),XX(~QC),100,'fill','s','MarkerFaceColor',col(1,:),'MarkerEdgeColor','k');
                hold off;
                ylabel('log(Peak Area)');
            else
                hold on;
                hqc = scatter(T(QC),XX(QC),100,'fill','o','MarkerFaceColor',col(2,:),'MarkerEdgeColor','k');
                hs = scatter(T(~QC),XX(~QC),100,'fill','s','MarkerFaceColor',col(1,:),'MarkerEdgeColor','k');
                hold off;
                ylabel('Peak Area');
            end;
            title(['Batch ',num2str(batch_num),'; Peak ',peak,'; m/z ',num2str(table2array(Peaks(peak,'Mass')),'%.5f'),' -- AFTER CORRECTION']);
            hold on;

            qcmean = nanmean(XX(QC));
            qcstd = nanstd(XX(QC));
            samplemean = nanmean(XX(~QC));
            samplestd = nanstd(XX(~QC));

            axis([min(T)-1, max(T)+1, min([nanmin(XX),samplemean-(2*samplestd)])*0.99, max([nanmax(XX),samplemean+(2*samplestd)])*1.01]);
            
            
            plot([0,r+1]',[samplemean+(2*samplestd),samplemean+(2*samplestd)],'--k');
            text(r+1,samplemean+(2*samplestd),'+2SD (Sample)');
            plot([0,r+1]',[samplemean-(2*samplestd),samplemean-(2*samplestd)],'--k');
            text(r+1,samplemean-(2*samplestd),'-2SD (Sample)');
            plot([0,r+1]',[qcmean+(2*qcstd),qcmean+(2*qcstd)],'--b');
            text(r+1,qcmean+(2*qcstd),'+2SD (QC)');
            plot([0,r+1]',[qcmean-(2*qcstd),qcmean-(2*qcstd)],'--b');
            text(r+1,qcmean-(2*qcstd),'-2SD (QC)');
            plot([0,r+1]',[qcmean,qcmean],'-.b');
            text(r+1,qcmean,'Mean QC');
            plot([0,r+1]',[samplemean,samplemean],'-.k');
            text(r+1,samplemean,'Mean Sample');
            
            hold off;
            
        end
        
        subplot(2,3,4)
        plot(params,R.cvMse,'b-x')
        xlabel('Param values');
        ylabel('cvMSE');
        axis tight;
        
        subplot(2,3,1)
        hold on;
        text(-0.5,1.0,['Peak: ',peak],'FontSize',10)
        text(-0.5,0.9,['MPA: ',num2str(R.MPA,'%8.0f')],'FontSize',10)
        text(-0.5,0.8,['QCwindow: ',num2str(R.QCwindow)],'FontSize',10)
        text(-0.5,0.7,['Optimal Param value: ',num2str(R.opt_param)],'FontSize',10)
        text(-0.5,0.6,['Leadin: ',num2str(R.leadin)],'FontSize',10)
        text(-0.5,0.5,['Leadout: ',num2str(R.leadout)],'FontSize',10)
        text(-0.5,0.4,['maxdist: ',num2str(R.maxdist)],'FontSize',10)
        text(-0.5,0.3,['MAXwindow: ',num2str(R.MAXwindow)],'FontSize',10)
        text(-0.5,0.2,['Total expected QCs: ',num2str(sum(QC))],'FontSize',10)
        text(-0.5,0.1,['Total actual QCs: ',num2str(R.numQCs)],'FontSize',10)
        if R.Poly == 1
            text(-0.5,-0.0,'Correction method: Poly');
        elseif R.Poly == 0
            text(-0.5,-0.0,'Correction method: Linear');
        elseif R.Poly == -2
            text(-0.5,-0.0,'Correction method: Mean ');
        else
            text(-0.5,-0.0,'Correction method: NONE ');
        end
        text(-0.5,-0.1,['~%RSD [BEFORE]AFTER: [',num2str(BRSD,'%3.1f'),'] ',num2str(R.RSD,'%3.1f')],'FontSize',10)
        text(-0.5,-0.2,['SNR(db) [BEFORE]AFTER : [',num2str(BSNR,'%3.1f'),'] ',num2str(R.SNR,'%3.1f')],'FontSize',10)
        axis off
        hold off;
        
        set(gcf, 'PaperType', 'A3');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0.25 2.07 10.50 4.37]);
        mz_save = num2str(table2array(Peaks(peak,'Mass')),'%.5f');
        mz_save = strrep(mz_save,'.','_');
        jpeg_out = [outdir,filesep,peak,'_','Batch_',num2str(batch_num),'_','_mz_',mz_save,'_QCRSC.jpg'];
        print(gcf,'-djpeg','-r300',jpeg_out);
        if ~isempty(jpeg_out)
            fns_jpeg{fign} = jpeg_out;
        end
    end
end
