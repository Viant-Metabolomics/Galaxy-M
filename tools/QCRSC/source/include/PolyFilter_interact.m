function [] = PolyFilter_interact(config_filename)

warning off;

[C] = load_config(config_filename);
[Meta,Data,Peaks] = load_data(C);
list = Peaks.Properties.RowNames;


Batchlist = unique(Meta.Batch);

display('');
display('-----------PolyFilter Interact Version 1.0 ----------------');
display('');

while 1
    inval = input('Which Batch would you like to investigate?','s');
    if strcmpi(inval,'end'), return; end;
    if strcmpi(inval,'q'), return; end;
    if strcmpi(inval,'quit'), return; end;
    inval = str2double(inval);
    if ~isempty(intersect(Batchlist,inval))
        batch_num = inval;
        break;
    else
        display('Batch number does not exist');
    end;
end

while 1

inval = input('Which Peak would you like to filter?','s');

if strcmpi(inval,'end'), return; end;
if strcmpi(inval,'q'), return; end;
if strcmpi(inval,'quit'), return; end;

if strcmpi(inval,'B')
    while 1
        inval = input('Which Batch would you like to investigate?','s');
        if strcmpi(inval,'end'), return; end;
        inval = str2double(inval);
        if ~isempty(intersect(Batchlist,inval))
            batch_num = inval;
            break;
        else
            display('Batch number does not exist');
        end;
    end
    continue;
end
    

peak = inval;
    
if strcmpi(peak,'R')
    pp = randperm(numel(list));
    peak = list{pp(1)};
end

index = ismember(list,peak);

if sum(index) == 0
    display(['Fatal Error: peak ',peak,' does not exist']);
    continue;
end;

if sum(index) > 1
    display(['Fatal Error: more than one peak called',peak]);
    continue;
end;

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
try
mdl = fit(Tsample,Xsample,C.polyFunc,opts);
catch
    display('Could not fit a curve');
    continue;
end

col = [255 10 71; 51 204 255]./255;
Label1 = cell(size(QC));
for i = 1:numel(QC)
    if QC(i) == 1, Label1{i} = 'QC'; else Label1{i} = 'Sample'; end
end

scrsz = get(0,'ScreenSize');
figure('Position',[(scrsz(3)*1/40) scrsz(4)/8 scrsz(3)*0.6 scrsz(4)*0.4]);
hold on;

CI = predint(mdl,T,C.polyCI,'observation','off');


hqc = scatter(T(QC),X(QC),100,col(2,:),'o','fill','MarkerEdgeColor','k');
hs = scatter(T(~QC),X(~QC),100,col(1,:),'s','fill','MarkerEdgeColor','k');

h = plot(T,CI,'m--');

%h = plot(mdl,'predobs',C.polyCI);


set(h,'LineWidth',1.5);
if strcmp(C.peakTransform,'log')
   xlabel('Injection Order \rightarrow','FontSize',14);
   ylabel('log(Peak Area)','FontSize',14);
else
   xlabel('Injection Order \rightarrow','FontSize',14);
   ylabel('Peak Area','FontSize',14);
end;


grid on;
title(['Batch ',num2str(batch_num),'; Peak ',peak],'FontSize',14);
xlim([min(T)-1, max(T)+1]);

%CI = predint(mdl,T,C.polyCI);


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

[pathstr,~,~] = fileparts(C.DataFile);

if C.plotflag == 1
    set(gcf, 'PaperType', 'A3');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 2.07 10.50 4.37]);
    print(gcf,'-djpeg','-r300',[pathstr,filesep,C.ProjectName,'_',peak,'_B',num2str(batch_num),'_polyfit_',num2str(C.polyFunc),'_',num2str(C.polyCI),'.jpg']);
elseif C.plotflag == 0
    reply = input('Do you want yo save this figure? Y/N [N]: ', 's');
    if strcmpi(reply,'Y')
        set(gcf, 'PaperType', 'A3');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0.25 2.07 10.50 4.37]);
        print(gcf,'-djpeg','-r300',[pathstr,filesep,C.ProjectName,'_',peak,'_B',num2str(batch_num),'_polyfit_',num2str(C.polyFunc),'_',num2str(C.polyCI),'.jpg']);
    end
else
    set(gcf, 'PaperType', 'A4');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 1.80 9.5 4.66]);
    print(gcf,'-dpsc2','-r300','-append',[pathstr,filesep,C.ProjectName,'_',peak,'_B',num2str(batch_num),'_polyfit_',num2str(C.polyFunc),'_',num2str(C.polyCI),'.ps']);
end

end

warning on;

end

