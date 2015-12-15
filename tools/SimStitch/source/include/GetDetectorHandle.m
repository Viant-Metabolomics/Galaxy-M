function [hDetector] = GetDetectorHandle(hRaw)
%gets the detector handle for the raw handle

%constants
MAX_DETECTORTYPES = 5;

%list all possible Detector Types
hRawInfo = hRaw.get('RawInfo');
DetectorCount = hRawInfo.get('DetectorCount');
if DetectorCount > 1
    disp(['Number of detectors in .raw: ', num2str(DetectorCount)]);
    %DetectorsOfTypeCount
    %there are 6 possibilities to detectortype. seems
    %XNo-Device is not an option
    for i = 0:(MAX_DETECTORTYPES-1)
        DetectorsOfTypeCount(i+1) = hRawInfo.get('DetectorsOfTypeCount',i);
    end
    disp('Detectors Of Type Count: '); disp(DetectorsOfTypeCount);
    detectorTypes = hRaw.set('Detector');
    disp('Detector Types: '); disp(detectorTypes);
    reply = input('Which Detector Type (0 to 5)? [0]:');
    if isempty(reply)
        reply = 0;
    end
else
    reply = 0;
end

% XDetectorRead based on detector type choice
try
    hDetector = hRaw.get('Detector',reply,1);
catch
    Xerrors(lasterr);
    error('Error in GetDetectorHandle');
end

