function [peaksy, peaksh, peaksf, peaksr] = KCeInterpolate(y,f,e,thresh, P, Instrument,noise,base,noise_input)
% INPUT
% y, f:     spectral data points
% e:        KC exponent
% thresh:   min peak height to return
% 
% OUTPUT
% see above

% find indices of maxima
%  limitations:
%  maxima at postion (1) or (end) are not included
%  exact adjacent data counted as 2 peaks eg [0 1 1 1 0] -> peaks at [2 4]
peaksi = maxima(y,length(y)) + 1;
N = length(peaksi);

peaksy = zeros(1,N);
peaksh = zeros(1,N);
peaksf = zeros(1,N);
peaksr = zeros(1,N);

% loop through each peak
fprintf('Progress:  0%%');
tic;
n = 0;
for i=1:N
    %fit: 3 data points
    j = (peaksi(i)-1):(peaksi(i)+1);
    X = [];
    dfoff = f(j(2));   %avoid rounding errors in inversion by using relative values
    foff = f(j) - dfoff;
    yy = y(j);
    X(:,1) = (foff.^2).';
    X(:,2) = foff.';
    X(:,3) = ones(3,1);
    Y = (yy.^(1/e)).';
%         Z = X\Y    %Gaussian elimination
    [U,S,V] = svd(X);
    Z = (V*inv(S)*U')*Y;
    a = Z(1);
    b = Z(2);
    c = Z(3);
    %fitted peak centre
    x0 = foff(2) - abs(foff(2)-foff(1))/2 * (yy(3)-yy(1)) / (yy(1) - 2*yy(2) + yy(3));
    %find peak height
    y0 = (a*x0.^2 + b*x0 + c)^e;
    
    
    % estimate peak SNR
    calculate = 0;
    if thresh %Only filter when thresh > 0? 
        if noise_input.ESTIMATE == 1
            % Test is peak height is above noise threshold
            if y0 > thresh
                calculate = 1;
            end
        else
            % Test if peak SNR is above threshold
            Peak_Noise = interp1(f,noise,x0+dfoff);
            Peak_Base = interp1(f,base,x0+dfoff);
            SNR_value = (y0 - Peak_Base)./Peak_Noise;
            if SNR_value > noise_input.MIN_SNR
                calculate = 1;
            end
        end
    else
        calculate = 1;
    end
    
    % is peak wanted?
    if calculate
        % yes
        n = n+1;
        peaksf(n) = x0 + dfoff;
        peaksh(n) = y0;
        % calculate the resolution FWHM
        q = sqrt(b^2-4*a*(c-exp(log(y0/2)/e)));
        peaksr(n) = f2mz(x0+dfoff,P, Instrument, []) / abs(f2mz((-b+q)/(2*a)+dfoff,P, Instrument, []) - f2mz((-b-q)/(2*a)+dfoff,P, Instrument, [])); 
        % accurately integrate using KCe interpolation
        q = 4*a*c - b^2;
        qq = sqrt(-q);
        zeroX(1) = (-b+qq)/(2*a);
        zeroX(2) = (-b-qq)/(2*a);
        if length(zeroX)~=2, disp(zeroX); error('Check zero-crossings'); end
        %sort limits
        if zeroX(1)>zeroX(2), zeroX = fliplr(zeroX); end
        %integrate using maple - slow
%         peaksd(n) = (double(int(ys,xs,zeroX(1)-foff,zeroX(2)-foff)));
        %integrate using indefinite solution - as per paper- fast
        z = e - 0.5;
        k = 4*a/q;          %NOTE MISPRINT IN PAPER!!!! ARRRGGHGHHH!!!  Took me hours to find that!
        ya = XRootXIntegral(z,k,a,b,c,zeroX(2),q) - XRootXIntegral(z,k,a,b,c,zeroX(1),q);
        %may be small imaginary component due to rounding errors (allow 0.1%)
        if abs(imag(ya)/real(ya)) > 0.001, error('Check imag components of peak area'); end
        peaksy(n) = real(ya);
    end

    % progress
    if floor(toc) || (i==N) % update every second
        fprintf('\b\b\b');
        fprintf('%2.0f%%',floor(i/N*100));
        tic;
    end
end
fprintf('\n');
% resize matrices
peaksy = peaksy(1:n);
peaksh = peaksh(1:n);
peaksf = peaksf(1:n);
peaksr = peaksr(1:n);
end

