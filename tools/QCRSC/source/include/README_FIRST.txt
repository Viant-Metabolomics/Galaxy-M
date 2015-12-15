%% 4 sections

% First examples of how to use the interactive scripts. Always use these first to get an idea of how the configuration settings perform. Then the full (and possibly slow) correction algorithm once you have optimized the config file. Usually the standard config file is OK

% The config.txt file is the most important file. Hopefully the first section is obvious.
% The sections regarding filter and correction and cleaning are a bit more complex. 
% Particularly tricky is the last section. Read our paper, and play around with values.


% 1. Polyfilter INTERACTIVE SCRIPT (prescreening obvious outliers)
% 2. QCRSC INTERACTIVE SCRIPT (looking at the correction for peak N of batch M)
% 3. Before After INTERACTIVE SCRIPT (looking at before after plots for multiple batches - % good for checking that it works.
% 4. Full QCRSC SCRIPT (multiple report files produced - one of the most the important % files is the “peaks removed” file - always check this)/ 

%%%%%%%%%%%%%% Polyfilter INTERACTIVE SCRIPT %%%%%%%%%%%%%%%%%%
%%% r = Random peak
%%% b = change batch
%%% q = quit
%%% end = quit
%%% quit = quit


PolyFilter_interact('configx.txt');
Loading Configuration file ...
 
ProjectName = Test
DataFile = Data.csv
PeakFile = Peaks.csv
PolyFilter = 1
QCRSC = 1
DataClean = 1
report = 1
parallel = 1
saveopt = 1
zeroflag = 1
peakTransform = log
plotflag = 1
polyFunc = poly3
polyCI = 0.99
polyAction = ignore
operator = subtract
maxWindow = -1
searchRange = 0:0.5:7
MPAtol = 4
kill = 1
pdiff = 0.0001
QCdist = 3
MaxRSD = 25
MinSNR = 0
ployAction = ignore
KsiTol = 4
 
... done.
 
Loading DataFile: Data.csv
TOTAL SAMPLES: 35 TOTAL PEAKS: 1929 NUMBER OF BATCHES: 1
Loading PeakFile: Peaks.csv
Done!
     ''

-----------PolyFilter Interact Version 1.0 ----------------
     ''

Which Batch would you like to investigate?1
Which Peak would you like to filter?r
Which Peak would you like to filter?M20
Which Peak would you like to filter?q



%%%%%%%%%%%%%% QCRSC INTERACTIVE SCRIPT %%%%%%%%%%%%%%%%%%
%%% r = Random peak
%%% b = change batch
%%% q = quit
%%% end = quit
%%% quit = quit


QCRSC_interact('configx.txt');
Loading Configuration file ...
 
ProjectName = Test
DataFile = Data.csv
PeakFile = Peaks.csv
PolyFilter = 1
QCRSC = 1
DataClean = 1
report = 1
parallel = 1
saveopt = 1
zeroflag = 1
peakTransform = log
plotflag = 1
polyFunc = poly3
polyCI = 0.99
polyAction = ignore
operator = subtract
maxWindow = -1
searchRange = 0:0.5:7
MPAtol = 4
kill = 1
pdiff = 0.0001
QCdist = 3
MaxRSD = 25
MinSNR = 0
ployAction = ignore
KsiTol = 4
 
... done.
 
Loading DataFile: Data.csv
TOTAL SAMPLES: 35 TOTAL PEAKS: 1929 NUMBER OF BATCHES: 1
Loading PeakFile: Peaks.csv
Done!
     ''

-----------QCRSC Interact Version 1.0 ----------------
     ''

Which Batch would you like to investigate?1
Which Peak would you like to correct?r
 
Peak M609
MPA is 107000
QCwindow is 2.0909
Optimal Param value is 0.5
Leadin is 0
Leadout is 0
maxdist is 4
MAXwindow is 0
Correction method: Poly
Total expected QCs is 11
Total actual QCs is 11
%RSD is 3.6007
SNR(dB) is 13.7455
Which Peak would you like to correct?M10
 
Peak M10
MPA is 18049.3767
QCwindow is 2.0909
Optimal Param value is 2
Leadin is 0
Leadout is 0
maxdist is 7
MAXwindow is 2
Correction method: Poly
Total expected QCs is 11
Total actual QCs is 10
%RSD is 13.0524
SNR(dB) is 7.4862
Which Peak would you like to correct?B
Which Batch would you like to investigate?1
Which Peak would you like to correct?r
 
Peak M125
MPA is 79200
QCwindow is 2.0909
Optimal Param value is 2.5
Leadin is 0
Leadout is 0
maxdist is 4
MAXwindow is 0
Correction method: Poly
Total expected QCs is 11
Total actual QCs is 11
%RSD is 9.6798
SNR(dB) is 3.5884
Which Peak would you like to correct?quit



%%%%%%%%%%%%%% BEFORE & AFTER INTERACTIVE SCRIPT %%%%%%%%%%%%%%%%%%
%%% r = Random peak
%%% q = quit
%%% end = quit
%%% quit = quit


 
BeforeAfter_interact('configx.txt');
Loading Configuration file ...
 
ProjectName = Test
DataFile = Data.csv
PeakFile = Peaks.csv
PolyFilter = 1
QCRSC = 1
DataClean = 1
report = 1
parallel = 1
saveopt = 1
zeroflag = 1
peakTransform = log
plotflag = 1
polyFunc = poly3
polyCI = 0.99
polyAction = ignore
operator = subtract
maxWindow = -1
searchRange = 0:0.5:7
MPAtol = 4
kill = 1
pdiff = 0.0001
QCdist = 3
MaxRSD = 25
MinSNR = 0
ployAction = ignore
KsiTol = 4
 
... done.
 
Loading DataFile: Data.csv
TOTAL SAMPLES: 35 TOTAL PEAKS: 1929 NUMBER OF BATCHES: 1
Loading PeakFile: Peaks.csv
Done!
Loading DataFile: QCRSC_Test_Data.csv
TOTAL SAMPLES: 35 TOTAL Test_PEAKS: 1929 NUMBER OF BATCHES: 1
Loading PeakFile: QCRSC_Peak.csv
Done!
Which Peak would you like to plot?r
Which Peak would you like to plot?r
Which Peak would you like to plot?M100
Which Peak would you like to plot?quit

%%%%%%%%%%%%%% FULL QCRSC SCRIPT %%%%%%%%%%%%%%%%%%


[Data,Peaks,C] = QCRSC_X3('configx.txt');
Loading Configuration file ...
 
ProjectName = Test
DataFile = Data.csv
PeakFile = Peaks.csv
PolyFilter = 1
QCRSC = 1
DataClean = 1
report = 1
parallel = 1
saveopt = 1
zeroflag = 1
peakTransform = log
plotflag = 1
polyFunc = poly3
polyCI = 0.99
polyAction = ignore
operator = subtract
maxWindow = -1
searchRange = 0:0.5:7
MPAtol = 4
kill = 1
pdiff = 0.0001
QCdist = 3
MaxRSD = 25
MinSNR = 0
ployAction = ignore
KsiTol = 4
 
... done.
 
Loading DataFile: Data.csv
TOTAL SAMPLES: 35 TOTAL PEAKS: 1929 NUMBER OF BATCHES: 1
Loading PeakFile: Peaks.csv
Done!
 
-------- Performing PolyFilter --------
 
Processing Batch 1
... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... removing QC outliers
Saving data as: Filt_Test_Peak.csv & Filt_Test_Data.csv
 
-------- Performing QCRSC --------
 
Processing Batch 1
max_window 5
Correcting peaks ... this make take a long time\n
reporting is disabled due to parallel processing ... be very patient
Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.

....................
Saving data as: QCRSC_Test_Peak.csv & QCRSC_Test_Data.csv
 
-------- Performing Data Clean --------
 
266 Peaks Removed
Saving data as: Cleaned_Test_Peak.csv & Cleaned_Test_Data.csv



