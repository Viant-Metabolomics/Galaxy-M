function [] = QCRSC_interact(config_filename)

[C] = load_config(config_filename);


[pathstr,~,~] = fileparts(C.DataFile);

filename1 = [pathstr,filesep,'Filt_',C.ProjectName,'_data.csv'];

if exist(filename1,'file')
    while 1
    inval = input('Do you want to correct the PolyFiltered Data?','s');
    if strcmpi(inval,'end'), return; end;
    if strcmpi(inval,'y') || strcmpi(inval,'yes')
        display('Using Filterred data');
        C.DataFile = [pathstr,filesep,'Filt_',C.ProjectName,'_data.csv'];
        C.PeakFile = [pathstr,filesep,'Filt_',C.ProjectName,'_peak.csv'];
        break;
    else
        display('Using RAW data');
        break;
    end;
    end
end


[Meta,Data,Peaks] = load_data(C);
list = Peaks.Properties.RowNames;

Batchlist = unique(Meta.Batch);

display('');
display('-----------QCRSC Interact Version 1.0 ----------------');
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

inval = input('Which Peak would you like to correct?','s');

if strcmpi(inval,'end'), return; end;
if strcmpi(inval,'q'), return; end;
if strcmpi(inval,'quit'), return; end;

if strcmpi(inval,'B') || strcmpi(inval,'Batch')
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
    
if strcmpi(peak,'R') || strcmpi(peak,'Rand')
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
R.D_RATIO = mad(XXsample,1)/mad(XXqc,1);
BRSD = 100*1.4826*mad(Xqc,1)/nanmedian(Xqc);
BD_RATIO = mad(Xsample,1)/mad(Xqc,1);


if R.numQCs < 5
    R.RSD = 100*nanstd(XXqc)/nanmean(XXqc);
    R.D_RATIO = nanstd(XXsample)/nanstd(XXqc);
    BRSD = 100*nanstd(Xqc)/nanmean(Xqc);
    BD_RATIO = nanstd(Xsample)/nanstd(Xqc);
end






if C.report
    display(' ');
    display(['Peak ',peak]);
    display(['MPA is ',num2str(R.MPA)]);
    display(['QCwindow is ',num2str(R.QCwindow)]);
    display(['Optimal Param value is ',num2str(R.opt_param)]);
    display(['Leadin is ',num2str(R.leadin)]);
    display(['Leadout is ',num2str(R.leadout)]);
    display(['maxdist is ',num2str(R.maxdist)]);
    display(['MAXwindow is ',num2str(R.MAXwindow)]);
    if R.Poly == 1
        display('Correction method: Cubic');
    elseif R.Poly == 0
        display('Correction method: Linear');
    elseif R.Poly == -2
        display('Correction method: Mean');
    else
        display('Correction method: NONE ');
    end
    display(['Total expected QCs is ',num2str(sum(QC))]);
    display(['Total actual QCs is ',num2str(R.numQCs)]);
    display(['%RSD is ',num2str(R.RSD)]);
    display(['D RATIO is ',num2str(R.D_RATIO)]);
    
end

%---------------------------------------------------------------
col = [255 10 71; 51 204 255]./255;
Label1 = cell(size(QC));
for i = 1:numel(QC)
    if QC(i) == 1, Label1{i} = 'QC'; else Label1{i} = 'Sample'; end
end


scrsz = get(0,'ScreenSize');
figure('Position',[(scrsz(3)*1/40) scrsz(4)/8 scrsz(3)*0.6 scrsz(4)*0.4]);

% scrsz = get(0,'ScreenSize');
% figure('Position',[(scrsz(3)*1/40) scrsz(4)/6 scrsz(3)*6/10 scrsz(4)*6/10]);
subplot(2,3,2:3)
if strcmp(C.peakTransform,'log')
   h = gscatter(T,X,Label1,'k','so',10,'off','','log(Peak Area)');
   set(h(1), 'MarkerFaceColor',col(1,:));
   set(h(2), 'MarkerFaceColor',col(2,:));
else
   h = gscatter(T,X,Label1,'k','so',10,'off','','Peak Area');
   set(h(1), 'MarkerFaceColor',col(1,:));
   set(h(2), 'MarkerFaceColor',col(2,:));
end;



title(['Batch ',num2str(batch_num),'; Peak ',peak,' -- BEFORE CORRECTION']);
hold on;

if nansum(X) ~= 0

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
plot([0,r+1]',[qcmean,qcmean],'--g');
plot([0,r+1]',[samplemean,samplemean],'--g');


plot(T,z,'k^')
hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nansum(XX) ~= 0


subplot(2,3,5:6)

if strcmp(C.peakTransform,'log')
    XX = log(XX);
    h = gscatter(T,XX,Label1,'k','so',10,'off','Injection Order \rightarrow','log(Peak Area)');
    set(h(1), 'MarkerFaceColor',col(1,:));
    set(h(2), 'MarkerFaceColor',col(2,:));
else
    h = gscatter(T,XX,Label1,'k','so',10,'off','Injection Order \rightarrow','Peak Area');
    set(h(1), 'MarkerFaceColor',col(1,:));
    set(h(2), 'MarkerFaceColor',col(2,:));
end;
title(['Batch ',num2str(batch_num),'; Peak ',peak,' -- AFTER CORRECTION']);
hold on;


%r = max(T);


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
plot([0,r+1]',[qcmean,qcmean],'--g');
plot([0,r+1]',[samplemean,samplemean],'--g');

hold off;

end;


subplot(2,3,4)
plot(params,R.cvMse,'b-x')
xlabel('Param values');
ylabel('cvMSE');
axis tight;

subplot(2,3,1)
hold on;
text(0,1.0,['Peak: ',peak],'FontSize',10)
text(0,0.9,['MPA: ',num2str(R.MPA,'%8.0f')],'FontSize',10)
text(0,0.8,['QCwindow: ',num2str(R.QCwindow)],'FontSize',10)
text(0,0.7,['Optimal Param value: ',num2str(R.opt_param)],'FontSize',10)
text(0,0.6,['Leadin: ',num2str(R.leadin)],'FontSize',10)
text(0,0.5,['Leadout: ',num2str(R.leadout)],'FontSize',10)
text(0,0.4,['maxdist: ',num2str(R.maxdist)],'FontSize',10)
text(0,0.3,['MAXwindow: ',num2str(R.MAXwindow)],'FontSize',10)
text(0,0.2,['Total expected QCs: ',num2str(sum(QC))],'FontSize',10)
text(0,0.1,['Total actual QCs: ',num2str(R.numQCs)],'FontSize',10)
if R.Poly == 1
    text(0,-0.0,'Correction method: Poly');
elseif R.Poly == 0
    text(0,-0.0,'Correction method: Linear');
elseif R.Poly == -2
    text(0,-0.0,'Correction method: Mean ');
else
    text(0,-0.0,'Correction method: NONE ');
end
text(0,-0.1,['~%RSD [BEFORE]AFTER: [',num2str(BRSD,'%3.1f'),'] ',num2str(R.RSD,'%3.1f')],'FontSize',10)
text(0,-0.2,['D RATIO [BEFORE]AFTER : [',num2str(BD_RATIO,'%3.1f'),'] ',num2str(R.D_RATIO,'%3.1f')],'FontSize',10)
axis off
hold off;

if C.plotflag == 1
    set(gcf, 'PaperType', 'A3');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 2.07 10.50 4.37]);
    print(gcf,'-djpeg','-r300',[pathstr,filesep,C.ProjectName,'_B',num2str(batch_num),'_',peak,'_QCRSC.jpg']);
elseif C.plotflag == 0
    reply = input('Do you want yo save this figure? Y/N [N]: ', 's');
    if strcmpi(reply,'Y')
        set(gcf, 'PaperType', 'A3');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0.25 2.07 10.50 4.37]);
        print(gcf,'-djpeg','-r300',[pathstr,filesep,C.ProjectName,'_B',num2str(batch_num),'_',peak,'_QCRSC.jpg']);
    end
else
    set(gcf, 'PaperType', 'A4');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 1.80 9.5 4.66]);
    print(gcf,'-dpsc2','-r300','-append',[pathstr,filesep,C.ProjectName,'.ps']);
end

end
end

