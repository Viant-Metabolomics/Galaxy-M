function [rms,fitSigma,fRange, numPoints] = GetNoiseLevel_0_5(data,f,exregions, P, Instrument)
% This function estimates the noise level of data, assuming the noise
% follows a Rayleigh distribution and that the data is predominantly noise.
% A section of test data is iteratively picked and fitted to the Rayleigh
% distribution.  The threshold (on the PDF) corresponding to the presence
% of a single data value above that threshold is determined, and this is
% compared to the actual maximum value in the test range selected.  If the
% actual maximum is higher, the region contains signal and the next region
% is selected.
% The test region is first moved up the data by a jump of delta, if no
% noise region is found and the end of the data is reached, the test region
% is reduced in size by rangeScale, and the process iterates.
%
%% Details ML fit

% fit_ML_rayleigh - Maximum Likelihood fit of the rayleigh distribution of i.i.d. samples!.
%                  Given the samples of a rayleigh distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(r)=r*exp(-r^2/(2*s))/s
%    with parameter: s
%
% format:   result = fit_ML_rayleigh( x,hAx )
%
% input:    x   - vector, samples with rayleigh distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      s   - fitted parameter
%                      CRB - Cram?r-Rao Bound for the estimator value
%                      RMS - RMS error of the estimation 
%                      type- 'ML'
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The rayleigh destribution is given by:
%
%    p(x,s) = x / s * exp(-x^2/(2*s))
%
%    where x are the samples which distribute by the function p(x,s)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    f(Xn,s)   = Xn / s * exp(-Xn^2/(2*s))
%    L(s)      = f(X,s) = product_by_n( f(Xn,s) )
%              = PI(Xn) * (s^(-N)) * exp( -sum(Xn^2)/(2*s) )
%    log(L(s)) = sum(log(Xn)) - N*log(s) - sum(Xn^2)/(2*s)
%     
%    The maximum likelihood point is found by the derivative of log(L(s)) with respect to "s":
%
%    diff(log(L(s))) = -N/s + sum(Xn^2)/(2*s^2) = (N/s^2) * ( sum(Xn^2)/(2*N) - s ) = 
%                    = J(s) * (s_estimation - s)  
%
%    Therefore, the (efficient) estimator is given by:
%
%               s = sum( Xn^2 ) / (2 * N)
%
%    The Cram?r-Rao Bound for this estimation is:
%
%               VAR( s ) = 1/J(s) = (s^2)/N
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

%% INPUTS
% data, f: mass spectrum data points
% exregions: frequency ranges to exclude from search ([range1;range2;...])

% OUTPUTS
% rms: rms noise level (of actual data)
% fitSigma: Rayleigh parameter (of fitted data)
% fRange: frequency range [from,to] over which noise is measured

% version 0.1, 20/Sep/07 - from GetNoiseThreshold_0_1, this function now simply returns the noise level
%                          (ie does not calculated an adjusted threshold - this moved to Stitch).
%              16/Jan/08 - displays m/z range over which noise calculated.
% version 0.2  06/Feb/08 - returns freq range over which noise calculated.
%              25/Feb/08 - minor changes, added parameters
% version 0.3  28/Feb/08 - catches case where noise region overlaps with
%                           region of "known" noise and searches for
%                           alternative
%              10/Mar/08 - changed delta to "delta = round(numPoints/10);"
%                           from "delta = round(numPoints/2);"
%
% 24/04/2015 (JE) combined GetNoiseLevel_0_5 and fit_ML_rayleigh code.
% (should be faster than loading fit_ML_rayleigh n times in the main while
% loop for GetNoiseLevel_0_5
%% Code
MAXPOINTS = 10000;   % maximum (starting) number of data points for parameter estimation
MINPOINTS = 50;    % minimum number of data points for parameter estimation

numPoints = min([MAXPOINTS length(data)]); % starting number of data points for noise estimation
delta = round(numPoints/100); % jump in index
rangeScale = 1.01;   % divisor for reducing the length of test range

%find a suitable range data with no real peaks
done = 0;   % flag
iStart = 1; % starting index
while ~done
    %select data points
    testData = data(iStart:iStart+numPoints-1);
    testf = f(iStart:iStart+numPoints-1);
    
    maxVal = max(testData); % Highest peak intensity in selection region; this intensity is later compared to the estimated maximum noise level
    %fit Rayleigh distribution (code replaceds ray = fit_ML_rayleigh(testData); )
    testData_x = real(testData(:));                 % should be column vectors !
    N       = length(testData_x);
    s       = sum(testData_x.^2)/(2*N);
    CRB     = (s^2)/N;
    [n,x_c] = hist( testData_x,100 );
    n       = n / sum(n*abs(x_c(2)-x_c(1)));
    y       = x_c.*exp(-x_c.^2/(2*s))/s;
    RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );
    ray = struct( 's',s,'CRB',CRB,'RMS',RMS,'type','ML' );
    
    sigma = sqrt(ray.s);
    maxLim = sigma*sqrt(-2 * log(1/numPoints));

    %move test region along by delta
    iStart = iStart + delta;
    
    %check if noise region includes exregions
    if find(min(exregions,[],2) < max(testf) & max(exregions,[],2) > min(testf))
        inEx = 1;
    else
        inEx = 0;
    end

    if (maxVal <= maxLim) && ~inEx
        % found region of noise
        fprintf('Found signal-free region over %d data points (%.2f-%.2fm/z), sigma = %0.2f\n',numPoints,f2mz(testf(end),P,Instrument),f2mz(testf(1),P,Instrument),sigma); %Update for other machines
        % debug output
        fprintf('Sigma from mean = %0.2f, rms error = %.2f\n',mean(testData)/sqrt(pi/2),ray.RMS);
        done = 1;
    elseif (iStart + numPoints > length(data))
        % got to end of data: decrease test region length and reset to
        % start of data
        iStart = 1;
        numPoints = round(numPoints/rangeScale);
        delta = round(numPoints/2);
        if numPoints < MINPOINTS
            plot(f2mz(f,P,instrument,[]),data); %Update for other machines
            error('Can''t find long enough region with no signal: check spectrum');
        end
    end    
end

rms = norm(testData)/sqrt(numPoints);
fitSigma = sigma;
fRange = [min(testf) max(testf)];
end

