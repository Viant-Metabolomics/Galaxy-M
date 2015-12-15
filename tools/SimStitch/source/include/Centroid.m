function [peaksy, peaksh, peaksf, peaksr] = Centroid(y,f,thresh,noise,base,noise_input)
%use the centroiding to return peaks
%uses only the maximum point and the two data points either side.

% indices of maxima
peaksi = maxima(y,length(y)) + 1;

% maximum data point for each peak
peaksm = y(peaksi);

% frequencies
peaksf=y(peaksi-1).*f(peaksi-1) + y(peaksi).*f(peaksi) + y(peaksi+1).*f(peaksi+1);
peaksf=peaksf./(y(peaksi-1) + y(peaksi) + y(peaksi+1));
%peaksf=data(peaksn-2).*f(peaksn-2) + data(peaksn-1).*f(peaksn-1) + data(peaksn).*f(peaksn) + data(peaksn+1).*f(peaksn+1) + data(peaksn+2).*f(peaksn+2);
%peaksf=peaksf./(data(peaksn-2) + data(peaksn-1) + data(peaksn) + data(peaksn+1) + data(peaksn+2));

% peak height (actually peak data maxima instead)
peaksh = peaksm;

% peak area (actually peak data maxima instead)
peaksy = peaksm;

% filter if required (include option NOT to filter)
if thresh
    if noise_input.estimate == 1
        idx = peaksh>thresh;
        peaksy = peaksy(idx);
        peaksf = peaksf(idx);
        peaksh = peaksh(idx);
    else
        
    end
end

% peak resolution: set to zero
peaksr = zeros(size(peaksy));

fprintf('\n');
end

