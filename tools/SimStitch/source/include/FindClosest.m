function [idxTarget,idxIn] = FindClosest_0_3(target,inputList,maxErr,errType)
%find THE (if present) closest target value to EACH inputList value
%ie. length(idxTarget) = length(idxIn) = length(inputList)
%returns the indices (if found) of the points in the inputList and the
%target list
%maxErr is maximum distance allowed from points
%errType = 'abs' or 'ppm' or 'off' for maxErr given as absolute or ppm or
%no maxErr.

%Updates:
%09/05/07, version 0.1: lists can be unsorted
%26/May/07, v0.1:   update for length(target)=1 and length(inputList)=1
%07/Dec/07, v0.2:   now returns only THE closest value, provided it is within
%                   the maxErr distance. Previously returned ALL closest values!
%18/Feb/08  v0.3:   minor bug correction in event only single target value
%                   (NB v0.2 seems to not have been renamed)

%check for empty input
if (isempty(inputList) | isempty(target)); warning('Empty input'); idxTarget=[]; idxIn=[]; return; end;

%sort inputList
[inputList,idxSortInput] = sort(inputList.');
inputList = inputList.';

%sort target
[target,idxSortTarget] = sort(target.');
target = target.';

%find index of closest in target list to each input
if length(target)==1
    ti = ones(size(inputList)); % each target is closest to input
else
    % interpolate inputList onto target
    ti = interp1(target, 1:length(target), inputList, 'nearest', 'extrap');
end

%find those with distance error <= maxErr
switch errType
    case 'off'
        ii = 1:length(inputList);
    case 'abs'
        e = abs(inputList - target(ti));
        ii = find(e <= maxErr);
        ti = ti(ii);
    case 'ppm'
        e = abs(inputList - target(ti)) ./ target(ti);
        ii = find(e <= maxErr/1e6);
        ti = ti(ii);
    otherwise
        error('Format choise not recognised');
end

%return indices wrt original (unsorted) input
idxTarget = idxSortTarget(ti).';
idxIn = idxSortInput(ii).';
end

