function vOutput = Apodize(vInput,method,QUIET)
% This function takes a vector input and outputs the same vector,
% multiplied by the apodization window function



%Determine the number of input data points
n = length(vInput);

%Initialize the vector
vFunc = linspace(0,n-1,n);

switch method
    case {0,'none'}
        if ~QUIET, disp('No apodisation'); end
        vOutput = vInput;
        return
    case {1,'hanning'}
        %Calculate the hanning funtion
        if ~QUIET, disp('Applying Hanning...'); end
        vFunc = .5*(1-cos(2*pi*vFunc/(n-1)));
    case {1,'hamming'}
        %Calculate the hamming funtion
        if ~QUIET, disp('Applying Hamming...'); end
        vFunc = 0.54 - 0.46*cos(2*pi*vFunc/(n-1));
    case {3,'blackman-harris'}
        if ~QUIET, disp('Applying Blackmann-Harris...'); end
        vFunc = 0.42323 - 0.49755*cos(2*pi*vFunc/(n-1)) + 0.07922*cos(4*pi*vFunc/(n-1));
    otherwise
        disp('ERROR: Choice not recognised! Default: no apodisation selected.');
        vOutput = vInput;
        return
end

%Output the result
vOutput = vInput.*vFunc;
if ~QUIET, disp('  ...done'); end
return
