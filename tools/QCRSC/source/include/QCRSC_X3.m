function [Data,peaks,Removed,C] = QCRSC_X3(config_filename)

[C] = load_config(config_filename);
Removed = [];

[Meta,Data,peaks] = load_data(C);

[pathstr,~,~] = fileparts(C.DataFile);

if C.PolyFilter
    if C.report
        display(' ');
        display('-------- Performing PolyFilter --------');
        display(' ');
    end
    [Data,peaks] = PolyFilter(Meta,Data,peaks,C);
    if C.saveopt
        if C.report, display(['Saving data as: Filt_',C.ProjectName,'_peaks.csv & Filt_',C.ProjectName,'_Data.csv']); end
        saveData(Meta,Data,peaks,[pathstr,filesep,'Filt_',C.ProjectName]);
    end
end

if C.QCRSC
    if C.report
        display(' ');
        display('-------- Performing QCRSC --------');
        display(' ');
    end
    [Data,peaks] = QCRSC_ALL(Meta,Data,peaks,C);
    if C.saveopt 
        if C.report, display(['Saving data as: QCRSC_',C.ProjectName,'_peaks.csv & QCRSC_',C.ProjectName,'_Data.csv']); end
        saveData(Meta,Data,peaks,[pathstr,filesep,'QCRSC_',C.ProjectName]);
    end
end

if C.DataClean
    if C.report
        display(' ');
        display('-------- Performing Data Clean --------');
        display(' ');
    end
    [Data,peaks,Removed] = DataClean(Meta,Data,peaks,C);
    if C.saveopt
        if C.report, display(['Saving data as: Cleaned_',C.ProjectName,'_peaks.csv & Cleaned_',C.ProjectName,'_Data.csv']); end
        saveData(Meta,Data,peaks,[pathstr,filesep,'Cleaned_',C.ProjectName]);
    end
end

if ~C.saveopt
    if C.report, display(['Saving data as: X3_',C.ProjectName,'_peaks.csv & X3_',C.ProjectName,'_Data.csv']); end
    saveData(Meta,Data,peaks,[pathstr,filesep,'X3_',C.ProjectName]);
end

end

