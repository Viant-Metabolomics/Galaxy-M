function [data,f] = ProcessTrans(transient,params,QUIET)

% MODIFICATIONS
% 
% 13/Mar/08     Updated mz2f to new function definition



PLOT_ON = 0;

%remove DC offset
if ~QUIET, disp('Removing DC offset...'); end
transient = transient - mean(transient);
if ~QUIET, disp('  ...done'); end
if PLOT_ON, plot(transient,'k'); pause; end   %f in MHz

%perform apodisation
transient = Apodize(transient,'hanning',QUIET);
if PLOT_ON, plot(transient,'k'); pause; end   %f in MHz

%convert to freq. domain
n = length(transient);
f1 = mz2f(params.mzEnd,params); %lower freq.
f2 = mz2f(params.mzStart,params); %upper freq.
[data,f] = Transform(transient,params.BW,'fft','magnitude',params.zFills,f1,f2,n,QUIET);
if PLOT_ON, plot(f./1e6,data); hold all; axis('auto'); pause; end   %f in MHz

%baseline correction
data = BaseLine(data,'off', QUIET);
return
