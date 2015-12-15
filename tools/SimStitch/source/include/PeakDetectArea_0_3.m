function [peaksy, peaksh, peaksf, peaksr] = PeakDetectArea_0_3(y,f,mode,thresh, P, instrument,noise,base,noise_input)
%find the peaks in data
% 
%output is
%-peaks .y (peak area)
%       .h (height)
%       .f (centre frequency)
%       .r (peak resolution)
% 
%inputs are:
%-f         frequency axis in Hz, increasing values
%-y         corresponding data points
%-mode      peak detect mode
%-thresh    minimum peak height - below this, peaks will not be returned

%CHANGES
% 22/Jun/07 Filter out peaks with <3 non-zero data points, area
%           determination not working for these (insufficient data for
%           determining peak area and center anyway)
% 01/Aug/07 (KCEinterp) include data points either side of maxima even if below noise
% 6/Sep/07  v0.1    Now also calculates the FWHM for peaks above the noise
%           floor and code optimised
% 6/Feb/08  v0.2    Added possibility to not apply peak threshold
%                   Returns maximum (existing) data point for peaks
% 14/Feb/08         thresh is now the maximum allowed data point (not peak
%                   height)
% 01/Apr/08         Updated centroid option. Removed QUIET. Removed
%                   redundant peak detection options
% 08/Apr/08 v0.3    Returns peak height instead of max. data point, and slightly optimised

switch mode
    case 'centroid'
        %centroid - computers in MS p35
        % quick: doesn't really work but gives an approximate answer
        [peaksy, peaksh, peaksf, peaksr] = Centroid(y,f,thresh,noise,base,noise_input);
    case 'KCEinterp'   %ref 56
        KCexp = 5.5;
        % other values of the KC exponent can be used for different line shapes - see ref 56
        [peaksy, peaksh, peaksf, peaksr] = KCeInterpolate(y,f,KCexp,thresh,P,instrument,noise,base,noise_input);
    otherwise
        error('Mode choice not recognised!');
end
end

