function [MZ, MZbar, Y, Ybar, N] = Clusters2Matrices_0_2(Ci, MZc, Yc, minPeaks, fOptions)
% Processes the output from PeakListCluster: expands the cluster labelling
% into matrices MZ, Y, N.
% Also filters the output according to minPeaks and fOptions.
% Sorts the result by MZbar
% 
% INPUT
% Ci - cell vector (length C) of cluster indices for each spectrum
% MZc - cell vector (length C) of spectrum m/z's
% Yc - cell vector (length C) of spectrum y values
% minPeaks - minimum number of peaks for n-filtering
% fOptions - mFilt 1 = remove CLUSTERS containing overlapping peaks =
%                      peaks from same spectrum occurring multiple times
%                  0 = use right-most peak in cluster)
%            nFilt 1 = filter peaks by presence in minimum number of
%                      spectra
%                  0 = no minimum number of spectra
%            remCl 1 = remove clusters filtered out by mFilt of nFilt
%                  0 = keep clusters but set N to 0 for all peaks caught
%                      by filter)
%             sort 1 = sort output by MZbar
%                  0 = don't
% 
% OUTPUT
% MZ - matrix (size C by P) of m/z values for each of P clusters
% MZbar - mean m/z values (length P) for each cluster
% Y - ", y values
% Ybar - ", mean y values (counting zeros as zeros)
% N - ", binary value indicating peak present (0) or absent (1)
% 
% VERSIONS
% 02/Apr/08 v0.1    Initial version
% 03/Apr/08 v0.2    Filtered clusters removed if remCl = 1, added sort as
%                   option

%% quick check
if (length(Ci) ~= length(MZc)) || (length(MZc) ~= length(Yc))
    error('Input not consistent');
end
numSpectra = length(MZc);

%% create cluster matrices
C = unique([Ci{:}]);    % index of cluster numbers to use
MZ = zeros(numSpectra,length(C));
Y = zeros(numSpectra,length(C));
N = zeros(numSpectra,length(C));
for i=1:numSpectra
    % m/z values
    MZ(i,Ci{i}) = MZc{i};
    % intensity values
    Y(i,Ci{i}) = Yc{i};
    % peak present
    N(i,Ci{i}) = 1;
end
% at this point, any overlapping peaks (multiple peaks present in same
% sample in same cluster) will be included but represented by the RIGHT
% most peak (consequence of assigning same value in MZ to multiple values
% in MZc above).
MZbar = MeanNonZero(MZ,1);
Ybar = mean(Y,1);

%% filter multiple peaks out
% if the same replicate has multiple peaks in a cluster, discard that
% cluster and make user aware until a more suitable way to
% differentiate between overlapping peaks can be found (also, currently,
% overlapping peaks will be poorly quantified)
if fOptions.mFilt
    edges = [C(1)-0.5 C+0.5];
    Cr = [];  % clusters to remove
    for i=1:numSpectra
        % count of each cluster
        cc = hist(Ci{i},edges,2);
        % select clusters to remove
        idx = cc>1;
        Cr = [Cr edges(idx)+0.5];
    end
    % found any?
    if ~isempty(Cr)
        warning('Found %d cluster(s) containing overlapping peaks',length(Cr));
        % yes - remove overlapping clusters?
        if fOptions.remCl
            % yes - remove unwanted clusters from index
            C = setdiff(C,Cr);
        else
            % no - just set cluster N's to zero
            N(:,Cr) = 0;
        end
    end
end

%% filter peaks occurring in minimum number of spectra
% filter on all clusters then combine result with any earlier filtering
if fOptions.nFilt
    % find peaks present in less than minimum number of replicates
    idx = sum(N,1)<minPeaks;
    Cr = find(idx);
    % any?
    if ~isempty(Cr)
        fprintf('Found %d clusters(s) containing insufficient peaks\n',length(Cr));
        % yes - remove cluster?
        if fOptions.remCl
            % yes - remove unwanted clusters from index
            C = setdiff(C,Cr);
        else
            % no - just set cluster N's to zero
            N(:,Cr) = 0;
        end
    end
end

%% apply filters to matrices
if fOptions.remCl
    % remove clusters if required
    MZ = MZ(:,C);
    MZbar = MZbar(C);
    Y = Y(:,C);
    Ybar = Ybar(:,C);
    N = N(:,C);
end

%% sort
if fOptions.sort
    % sort by increasing m/z
    [MZbar, idx] = sort(MZbar,'ascend');
    Y = Y(:,idx);
    Ybar = Ybar(:,idx);
    MZ = MZ(:,idx);
    N = N(:,idx);
end

return
end

