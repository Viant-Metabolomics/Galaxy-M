function [hRaw] = GetRawHandle(fileName)
%returns a new RAW handle given a filename

% MODIFICATION HISTORY
% 14/Mar/08         Added error code output.

try
    hRaw = actxserver('XRaw.XRaw');
    hRaw.Open(fileName);
catch
    errCode = lasterr;
    errCode = hex2dec(errCode(end-3:end));
    fprintf('Error loading file: %s\n',fileName);
    switch errCode
        case 10001
            error('File not found.');
        case 10002
            error('Internal error: wrong dll version? Try running reg.m');
        case 10006
            error('Cannot access information in the specified file type.');
        case 10007
            error('The file is already open.');
        case 10008
            error('No FileName supplied, or the name of the file was not set up in another method.');
        case 10009
            error('An invalid parameter value has been supplied.');
        case 10012
            error('The file is of an incorrect type.');
        case 10014
            error('There is no valid data available.');
        case 10017
            error('The object is not valid.');
        otherwise
            rethrow(lasterror);
    end
end
end

