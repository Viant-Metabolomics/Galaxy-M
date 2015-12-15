function outData = BaseLine (inData, mode, QUIET)
%Baseline correct


switch mode
    case {0,'off'}
        if ~QUIET, disp('No baseline correction'); end
        outData = inData;
        return
    %case {1, 'on'}
        %disp('Applying baseline correction...');
    otherwise
        disp('ERROR: Mode choice not recognised!');
        outData = inData;
        return
end

if ~QUIET, disp('  ...done'); end
return

