function [FSout, PLout] = SegStitch_1_1(FS, PL, P, marginToCut, Instrument)
%stitches multiple spectra together, assuming A and B are the same for all
%spectra, using spline interpolation to make the data points monotonic
%assumptions:
%noise has been removed from the spectra (ie set to zero)
%limitations:
%the frequency points must be evenly spaced
%the spacing in the first range will determine the spacing of the stitched spectrum

% INPUTS:
% FS: array of spectra data structs, one for each SIM window
% PL: array of spectra peak structs, one for each SIM window
% P: array of parameters structs, one for each SIM window
% marginToCut: edge effects in each SIM window to remove during stitching

% OUTPUTS:
% FSout: single spectra data struct (stitched)
% PLout: single spectra peak struct (stitched)

%version history:
%9/5/07, v0.2    redundant modes removed (leaving spline only)
%                removed references to specdf
%18/6/07         dead region now dependent upon window mzStart
%19/6/07         included test if no overlapping regions left, in which case peaks are concat to save time
%                added check to remove any data points beyond the spectral mzStart or mzEnd
%06/02/08, v0.3  replaced spec.peaksn with spec.peaksFlag and spec.pNoise
%                only cut spectra at minima below the data threshold
%22/02/08        minor error-causing bug fix
%26/02/08, v1.0  Removed params, numSpec and mode parameters, tidied code
% 13/Mar/08         Added output of by how much target cut point is moved
% 14/Apr/08         Minor change to remove old specParams variable line 78
% 14/Apr/08 v1.1    Added edge effect to remove from start of first SIM
%                   window and end of last SIM window

numWindows = length(P);

% catch case only 1 SIM window
if numWindows == 1
    FSout = FS;
    PLout = PL;
    
    % 19-10-2015 Extra code nature protocols paper 
    
    % Spectral data
    mz = f2mz(FSout(1).f,P(1),Instrument);
    idx = find(mz > (P(1).mzEnd - marginToCut.MZMAX) | mz < (P(1).mzStart + marginToCut.MZMIN)); % Remove null region from start and end of window 
    sFields = fieldnames(FSout);
    for i=1:length(sFields)
        FSout.(sFields{i})(idx) = []; % Update FS
    end
    mz(idx) = [];
    idx2 = mz > marginToCut.BOUND(end) | mz < marginToCut.BOUND(1); % Remove parts of window outside bounds
    for i=1:length(sFields)
        FSout.(sFields{i})(idx2) = []; % Update FS
    end
    
    % Spectral peaks
    mz = f2mz(PLout(1).f,P(1),Instrument);
    idx = find(mz > (P(1).mzEnd - marginToCut.MZMAX) | mz < (P(1).mzStart + marginToCut.MZMIN));
    sFields = fieldnames(PLout);
    for i=1:length(sFields)
        PLout.(sFields{i})(idx) = [];
    end
    mz(idx) = [];
    idx2 = mz > marginToCut.BOUND(2) | mz < marginToCut.BOUND(1); % Remove parts of window outside bounds
    for i=1:length(sFields)
        PLout.(sFields{i})(idx2) = [];
    end

    return
end

% Remove windows outside bounds
remove = [];
for si = 1:numWindows
    if P(si).mzEnd < marginToCut.BOUND(1)
        remove = [remove, si];
    end
    
    if P(si).mzStart > marginToCut.BOUND(end)
        remove = [remove, si];
    end
end
FS(remove) = [];
PL(remove) =  [];
P(remove) = [];
numWindows = length(P);

% determine a common freq. for each spectrum, will form the freq points of the new stitched spectrum
% the deltaf should be taken from a CALIBRATED (ie freq. points have not been adjusted) spectrum
deltaf = 0;
for si = 1:numWindows
    if (P(si).ecald || P(si).icald)
        deltaf = FS(si).f(2) - FS(si).f(1);
        break;
    end
end
if deltaf == 0
    error('No calibrated SIM windows found');
end

%cut the "edge effect" from spectra
%select edge effect based on midpoint of overlap
sw_l = zeros(numWindows,1);
sw_h = zeros(numWindows,1);
sw_l(1) = marginToCut.MZMIN;
for si=2:numWindows
    mp = (P(si).mzStart + P(si-1).mzEnd)/2;
    for i=1:length(marginToCut.BOUND)-1
        if (marginToCut.BOUND(i) < mp && mp <= marginToCut.BOUND(i+1))
            sw_l(si) = marginToCut.START(i);
            sw_h(si-1) = marginToCut.END(i);
        end
    end
end
sw_h(end) = marginToCut.MZMAX;
for si=1:numWindows
    fprintf('Segment %d (%d-%dm/z): %.2f (lower edge), %.2f (higher edge)\n',...
        si,P(si).mzStart,P(si).mzEnd,sw_l(si),sw_h(si));
    %if after removing the edge effects the windows still overlap, unpreferentially increase the edge effects until no overlap exists
    %this is to avoid averaging spectral values
    if si<numWindows
        olap = P(si).mzEnd - P(si+1).mzStart - sw_h(si) - sw_l(si+1);
        if olap > 0
            fprintf('INCREASING AMOUNT OF REMOVED EDGE EFFECT to zero overlap:\n');
            sw_h(si) = sw_h(si) + olap /2;
            sw_l(si+1) = sw_l(si) + olap /2;
            fprintf('Segment %d (%d-%dm/z): %.2f (lower edge), %.2f (higher edge)\n',si,P(si).mzStart,P(si).mzEnd,sw_l(si),sw_h(si));
        elseif olap < 0
            error('Too much edge effect: missing data points');
        end
    end
    if P(si).dNthresh > 0
        %only cut at a minima in the data below the threshold (if noise
        %non-zero)
        zIdx = maxima(-1*FS(si).data,length(FS(si).data)) + 1;
        zIdx = zIdx(FS(si).data(zIdx) < P(si).dNthresh);
    elseif length(find(FS(si).data==0)) > 0
        % only cut at a zero (if data contains zero values)
        zIdx = find(FS(si).data==0);
    else
        % only cut at a minimum (if noise not set and no zeros)
        zIdx = maxima(-1*FS(si).data,length(FS(si).data)) + 1;
    end
    %higher frequency end
    fHiCut = mz2f(P(si).mzStart + sw_l(si),P(si),Instrument);
    idx = find(FS(si).f <= fHiCut);
    iH = max(intersect(idx, zIdx));
    % check haven't moved way beyond target cut point
    fprintf('Moved low cut point by %.2f m/z\n',...
        abs(f2mz(FS(si).f(max(idx)),P(si), Instrument)-f2mz(FS(si).f(max(iH)),P(si), Instrument)));
    %lower frequency end
    fLoCut = mz2f(P(si).mzEnd - sw_h(si),P(si),Instrument);
    idx = find(FS(si).f >= fLoCut);
    iL = min(intersect(idx, zIdx));
    % check haven't moved way beyond target cut point
    fprintf('Moved high cut point by %.2f m/z\n',...
        abs(f2mz(FS(si).f(min(idx)),P(si),Instrument)-f2mz(FS(si).f(max(iL)),P(si),Instrument)));
    %cut
    sFields = fieldnames(FS);
    for i=1:length(sFields)
        FS(si).(sFields{i}) = FS(si).(sFields{i})(iL:iH);
    end

    %peaks - high freq
    iH = max(find(PL(si).f <= fHiCut));
    %lower frequency end
    iL = min(find(PL(si).f >= fLoCut));
    %cut
    sFields = fieldnames(PL);
    for i=1:length(sFields)
        PL(si).(sFields{i}) = PL(si).(sFields{i})(iL:iH);
    end
end

% trim spectra to mzStart and mzEnd if necessary
for si=1:numWindows
    % spectral data
    mz = f2mz(FS(si).f,P(si),Instrument);
    idx = find(P(si).mzEnd < mz | mz < P(si).mzStart);
    if ~isempty(idx)
        fprintf('Trimming spec %d (%.2f-%.2f): removing spec data from %.6f-%.6f\n',...
            si,P(si).mzStart,P(si).mzEnd,min(mz(idx)),max(mz(idx)));
        sFields = fieldnames(FS);
        for i=1:length(sFields)
            FS(si).(sFields{i})(idx) = [];
        end
    end
    % spectral peaks
    mz = f2mz(PL(si).f,P(si),Instrument);
    idx = find(P(si).mzEnd < mz | mz < P(si).mzStart);
    if ~isempty(idx)
        fprintf('Trimming spec %d (%.2f-%.2f): removing peaks from %.6f-%.6f\n',...
            si,P(si).mzStart,P(si).mzEnd,min(mz(idx)),max(mz(idx)));
        sFields = fieldnames(PL);
        for i=1:length(sFields)
            PL(si).(sFields{i})(idx) = [];
        end
    end
    
    % 19-10-2015 remove part window outside of bounds
    % Spectral data
    mz = f2mz(FS(si).f,P(si),Instrument);
    idx2 = mz > marginToCut.BOUND(end) | mz < marginToCut.BOUND(1); % Remove parts of window outside bounds
    if ~isempty(idx2)
        sFields = fieldnames(FS);
        for i=1:length(sFields)
            FS(si).(sFields{i})(idx2) = []; % Update FS
        end
    end
    
    % Spectral peaks
    mz = f2mz(PL(si).f,P(si),Instrument);
    idx2 = mz > marginToCut.BOUND(end) | mz < marginToCut.BOUND(1); % Remove parts of window outside bounds
    if ~isempty(idx2)
        sFields = fieldnames(PL);
        for i=1:length(sFields)
            PL(si).(sFields{i})(idx2) = [];
        end
    end

end

% ====================== 19-10-2015 =====================================
% Sf seems to be unused?
%next determine what Sf, the f values for the final spectrum, will be
%assuming the frequency values are already cut
%find the maximum and minimum values of Sf
SfMin = min(FS(1).f);
SfMax = max(FS(1).f);
for si = 1:numWindows
    fMax(si) = max(FS(si).f);
    fMin(si) = min(FS(si).f);
    SfMax = max([SfMax fMax(si)]);
    SfMin = min([SfMin fMin(si)]);
end

%create Sf
Sf = SfMax:-deltaf:SfMin;
%create max to min then flip to avoid a /0 error
%also makes sense to construct based on highest freq. (1st range).
Sf = fliplr(Sf);

% 5/12/06 - bug fix... if the first segment has been aligned, need to remove the offset here such that the
% first segment (because stitching is referenced from the 1st seg) aligns with original data points for
% KCe interpolation to work.
m = Sf(1) / deltaf;
segShift = (round(m)-m)*deltaf;
Sf = Sf + segShift;
% =======================================================================

%if no regions are now overlapping, simply concatenate the remaining
%lengths of windows (spectrum may not have monotonous frequency values for
%peak-picking)
if issorted([FS(end:-1:1).f])
    fprintf('No overlap remaining: concatenating windows...');
    % Frequency spectra
    sFields = fieldnames(FS);
    for i=1:length(sFields)
        FSout.(sFields{i}) = [FS(end:-1:1).(sFields{i})];
    end
    % Peak list
    sFields = fieldnames(PL);
    for i=1:length(sFields)
        PLout.(sFields{i}) = [PL(end:-1:1).(sFields{i})];
    end
    %check consistency
    if ~issorted(FSout.f), error('Inconsistency in stitching'); end
    if ~issorted(PLout.f), error('Inconsistency in stitching'); end
    fprintf('done\n');
else
    error('Code not capable of stitching overlapping windows');
end
end

