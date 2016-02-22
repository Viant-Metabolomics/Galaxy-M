function [fns_jpeg] = PlotPolyFilterData(config_filename, peaks_to_plot, outdir)

[C] = load_config(config_filename);
[Meta,Data,Peaks] = load_data(C);
list = Peaks.Properties.RowNames;

Batchlist = unique(Meta.Batch);

fns_jpeg = cell(1,length(peaks_to_plot));
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
        batch_num = Batchlist(j);
        batch = ismember(Meta.Batch,batch_num);
        
        X = Data.(peak);
        X = X(batch);
        T = Meta.Order(batch);
        QC = Meta.QC(batch);
        
        if (C.zeroflag == 1),  X(X == 0) = NaN; end;
        if strcmp(C.peakTransform,'log'),  X = log(X); end;
        
        temp = isnan(X);
        
        X(temp) = [];
        T(temp) = [];
        QC(temp) = [];
        
        Xsample = X(~QC);
        Tsample = T(~QC);
        
        opts = fitoptions(C.polyFunc, 'Robust', 'Bisquare');
        
        mdl = fit(Tsample,Xsample,C.polyFunc,opts);
        
        col = [255 10 71; 51 204 255]./255;
        Label1 = cell(size(QC));
        for i = 1:numel(QC)
            if QC(i) == 1, Label1{i} = 'QC'; else Label1{i} = 'Sample'; end
        end
        
        scrsz = get(0,'ScreenSize');
        h = figure('Position',[(scrsz(3)*1/40) scrsz(4)/8 scrsz(3)*0.6 scrsz(4)*0.4]);
        set(h,'Visible','off');
        hold on;
        hqc = scatter(T(QC),X(QC),100,'fill','o','MarkerFaceColor',col(2,:),'MarkerEdgeColor','k');
        hs = scatter(T(~QC),X(~QC),100,'fill','s','MarkerFaceColor',col(1,:),'MarkerEdgeColor','k');
        h = plot(mdl,'predobs',C.polyCI);
        set(h,'LineWidth',1.5);
        if strcmp(C.peakTransform,'log')
            xlabel('Injection Order \rightarrow','FontSize',14);
            ylabel('log(Peak Area)','FontSize',14);
        else
            xlabel('Injection Order \rightarrow','FontSize',14);
            ylabel('Peak Area','FontSize',14);
        end
        grid on;
        title(['Batch ',num2str(batch_num),'; Peak ',peak,'; m/z ',num2str(table2array(Peaks(peak,'Mass')),'%.5f')],'FontSize',14);
        xlim([min(T)-1, max(T)+1]);
        
        CI = predint(mdl,T,C.polyCI);
        
        Xqc = X(QC);
        Tqc = T(QC);
        CIqc = CI(QC,:);
        CIsample = CI(~QC,:);
        
        cut_up = Xqc > CIqc(:,2);
        cut_low = Xqc < CIqc(:,1);
        QCcut = or(cut_up,cut_low);
        
        Xqc_cut = Xqc(QCcut);
        Tqc_cut = Tqc(QCcut);
        
        cut_up = Xsample > CIsample(:,2);
        cut_low = Xsample < CIsample(:,1);
        Scut = or(cut_up,cut_low);
        
        Xsample_cut = Xsample(Scut);
        Tsample_cut = Tsample(Scut);
        
        h2 = plot(Tqc_cut,Xqc_cut,'ko','MarkerSize',14,'LineWidth',1);
        h3 = plot(Tsample_cut,Xsample_cut,'k^','MarkerSize',14,'LineWidth',1);
        
        
        header = hs;
        leg = {'Samples'};
        
        if sum(QC)
            header = [header;hqc];
            leg = [leg;{'QCs'}];
        end
        
        header = [header;h];
        leg = [leg;{C.polyFunc};{[num2str(100*C.polyCI),'% Lower bounds']};{[num2str(100*C.polyCI),'% Upper bounds']}];
        
        if sum(QCcut)
            header = [header;h2];
            leg = [leg;{['QC Outlier (n = ',num2str(sum(QCcut)),')']}];
        end
        
        if sum(Scut)
            header = [header;h3];
            leg = [leg;{['Sample Outlier (n = ',num2str(sum(Scut)),')']}];
        end
        
        legend(header,leg, 'Location', 'NorthEastOutside','FontSize',14);
        box on;
        hold off;
        
        set(gcf, 'PaperType', 'A3');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        %set(gcf, 'PaperPosition', [0.25 2.07 10.50 4.37]);
        set(gcf, 'PaperPosition', [0.0 0.0 12.00 8.00]);
        mz_save = num2str(table2array(Peaks(peak,'Mass')),'%.5f');
        mz_save = strrep(mz_save,'.','_');
        jpeg_out = [outdir,filesep,peak,'_','Batch_',num2str(batch_num),'_','_mz_',mz_save,'_polyfit_',num2str(C.polyFunc),'_',num2str(C.polyCI),'.jpg'];
        print(gcf,'-djpeg','-r300',jpeg_out);
        if ~isempty(jpeg_out)
            fns_jpeg{k} = jpeg_out;
        end
    end
end
close all force
end
