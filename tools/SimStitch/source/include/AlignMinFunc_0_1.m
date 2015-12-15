function resid = AlignMinFunc_0_1(rL,sL,rR,sR,params)
% this function is the minimisation function for the alignment of
% overlapping SIM windows. The function is the sum of squared errors
% between the frequency points

% based on model that s = m.(r.f) + c {+a.((r.f).^2)} (not implemented)

% alignment could be weighted to the intensity of the peaks (ie multiply
% the error by the data)

% INPUTS:
% rL,sL,rR,sR: reference and subject (ie subject) spectra, containing
%   fields .y and .f: left and right (as applicable).  Align only non-empty
%   spectra
% params: [c,m,a] offset (0 order), scaling (1st order) and 2nd order
%   (not implemented) parameters. m and a parameters are optional

% Change log

% 03/Mar/08 v0.1    Initial version

c = params(1);
if length(params) == 2
    m = params(2) / 1e6;    % to improve sensitivity of fminsearch
else
    m = 1;
end

% number of peaks using to align
if isempty(rL)
    nL = 0;
else
    nL = length(rL.f);
end
if isempty(rR)
    nR = 0;
else
    nR = length(rR.f);
end

% apply alignment
if nL
    sL.f = m*sL.f + c;
end
if nR
    sR.f = m*sR.f + c;
end

% sum of squares
resid = 0;
if nL
    resid = resid + sum( (sL.f-rL.f).^2 );
end
if nR
    resid = resid + sum( (sR.f-rR.f).^2 );
end
end

