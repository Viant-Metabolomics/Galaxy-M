function [fdata,f] = Transform(tdata, BW, transform, mode, zeroFill, f1, f2, n, QUIET)
%Perform frequency domain transformation with options



switch transform
    case {0, 'fft'}
        if isempty(zeroFill)
            disp('WARNING: No zero-fill specified, Default: 1 zero-fill selected');
            zeroFill = 1;
        end
        if ~QUIET, disp(['Applying FFT with ',num2str(zeroFill),' zero fill(s), BW=',num2str(BW),'...']); end
        tic;
        n = (2^zeroFill) * length(tdata);
        fdata = fft(tdata, n);     %fft'd data
        fdata = fdata(1:n/2+1);     %first half of spectrum required
        fs = 2*BW;
        df = fs/n;  %freq step
        f = (0:n/2)*df;    %frequency scale
        %mzEps = 0;      %add m/z 1 margin overlap
        %f1 = f1 - mz2f(mzEps,A,B);
        %f2 = f2 + mz2f(mzEps,A,B);
        %chop frequency
        firstIdx = min(find(f>f1));
        lastIdx = max(find(f<f2));
        f = f(firstIdx:lastIdx);
        fdata = fdata(firstIdx:lastIdx);
    case {1, 'chirp'}
        disp(['Applying CZT with BW=',num2str(BW),'...']);
        tic;
        fs = BW*2;
        w = exp(-j*2*pi*(f2-f1)/(n*fs));
        a = exp(j*2*pi*f1/fs);
        fdata = czt(tdata,n,w,a);
        f = ((0:length(fdata)-1)*(f2-f1)/length(fdata)) + f1;
    otherwise
        disp('ERROR: Transform choice not recognised!');
        return
end

if ~QUIET, disp([' ...done (elapsed time: ',num2str(toc),'s)']); end

switch mode
    case {0,'magnitude'}
        if ~QUIET, disp('Calculating magnitude-mode...'); end
        fdata = abs(fdata); %magnitude mode
    otherwise
        disp('ERROR: Mode choice not recognised!');
        return
end

if ~QUIET, disp('  ...done'); end
return

