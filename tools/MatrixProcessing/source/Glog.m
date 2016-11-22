function Glog(dso_xml, outfile)
% -------------------------------------------------------------------
% Performs glog optimisation on the input dataset. If a dataset object,
% only the QCs should be tagged as 'included' so that only they will be
% used for optimisation. The output dataset object dso_out will contained
% the glog transformed spectra for all samples, not just QCs. 
%  
% If a peak intensity matrix is used, only the QCs should be included. The
% output matrix, dso_out will only contain the transformed QCs so the
% remaining samples should be transformed using lambda and the scaling
% factor. This can be performed for a matrix of samples as follows:
% 
%   scaled_samples = samples.*scaling_factor;
%   scaled_samples = scaled_samples.^power_factor;
%   glogged_samples = log(scaled_samples+ sqrt(power(scaled_samples,2)+lambda))
%
% 
% Usage: 
% [dso_out, lambda, scaling_factor]=glogopt_rld(dso);
%
%   inputs: 
%           dso: can be a matrix or a dataset object
%                   NB: if a dataset object, make sure only QC samples are
%                   'included'
%
%   outputs:
%           dso_out: the glog transformed data (will transform all spectra
%           even if only QCs were 'included' in input.
%           lambda = optimised or standard glog parameter used
% 
%   Code adapted by R L Davidson, 23rd November 2011
% -------------------------------------------------------------------

%% CHECK FOR PLS TOOLBOX IN PATH
if ~isdeployed
	try
	    dso = dataset(rand(10,100));       %attempt to create a dataset object
	    props = properties(dso);           %request the properties of said object
	catch err
	    if strcmp(err.identifier,'MATLAB:UndefinedFunction') %if dataset function not available then neither PLS Toolbox or Stats Toolbox are installed.
		
		%fid = fopen(messagefile,'w');
		sprintf('Matlab does not recognise dataset function. Neither PLSToolbox or Stats toolbox installed. Please amend and try again.');
		%fclose(fid);
		return

	    end
	end





	if ~isempty(props) %PLS Toolbox datasets have no properties at initiation whereas Matlab datasets have 2.
	    
	    % if here, Statistics Toolbox version has been used. need to move stats toolbox below pls toolbox.
	    clear dso
	    clear classes   %need to remove the statistics toolbox dataset class
	    
	    original_path = path; %save original path
	    rem = original_path;
	    pls_path = '';
	    rem_path = '';
	    while true
		[str,rem] = strtok(rem,pathsep);
		if isempty(str)
		    break
		elseif strfind(str,'pls_toolbox') %covers all pls_toolbox entries
		    pls_path = [pls_path,str,pathsep];
		else
		    rem_path = [rem_path,str,pathsep];
		end
	    end
	    
	    if ~isempty(pls_path) %check for no PLS toolbox installed
		path(pls_path,rem_path); %put PLS at the top!
		rehash pathreset;
		rehash toolboxreset;
		
		dso = dataset(rand(10,100));
		props = properties(dso);
		if ~isempty(props) % if rehash has not worked, quit.
		    %fid = fopen(messagefile,'w');
		    sprintf('Cannot appropriately rejig path. Please manually place PLSToolbox above Stats Toolbox in path.');
		    %fclose(fid);
		    path(original_path);
		    return
		end
		
	    else    % If no stats toolbox entries have been found in path there is a more serious problem.
		%fid = fopen(messagefile,'w');
		sprintf('PLS Toolbox not on path. Please Install and try again.');
		%fclose(fid);
		path(original_path);
		return
	    end
	    
	else
	    sprintf('PLS Toolbox dataset objects are available. Continuing.')
	    original_path = path;
	end
end


%% LOAD PLS DATASET OBJECT

[dso, name, source] = autoimport(dso_xml, 'xml');

%% DATA CHECKING
% Automatically include QC samples
[n,~] = size(dso);
QC_samples = [];
for i = 1:n
    if strcmp(dso.classid{1}(i),'QC') || strcmp(dso.classid{1}(i),'qc') || strcmp(dso.classid{1}(i),'Qc') || strcmp(dso.classid{1}(i),'qC')
        QC_samples = [QC_samples i];
    end
end
if isempty(QC_samples);
    warning('No QC samples detected, glog transformation should not be used');
else
    dso.include{1} = QC_samples;
end

%check for dataset object
if isobject(dso)
    x = dso.data;
    %test for those that are 'included'
    inc = dso.include{1};
    x = x(inc,:);
    x = x'; %for the purpose of the current implementation of the glog algorithm, the order of columns and rows has to be transposed. 
else
    x = dso'; %if the user passes a simple peak matrix matlab variable, they should be told to submit it in the rows=samples, cols=variables formation that is used by the PLS Toolbox dataset object. It then makes sense to transpose it here. Much easier to treat both the same (and better still to alter the GLOG algorithm below)
end

% open  'outfile' for status reporting
%fid = fopen(messagefile,'a');

%% PERFORM GLOG OPTIMISATION


% y0 is not currently used but can be used as an offset (I think).
y0=0;
N=size(x,2);
L=max(size(x));

%scale the data - required so that there is enough variance for the
%algorithm to optimise with. 
test = min(var(x'));
test2 = max(var(x'))./min(var(x'));
count = 0;
scaling_factor = 1;
power_factor = 1;
x_tmp = x;

while test2 <=1e10 || test <=1e5 
    if test<=1e5
        scaling_factor = scaling_factor*10;
    end
    if test2 <=1e10
        power_factor = power_factor+0.1;
    end
    x_tmp = x.*scaling_factor;
    x_tmp = x_tmp.^power_factor;
    count = count+1;
    test = min(var(x_tmp'));
    test2 = max(var(x_tmp'))./min(var(x_tmp'));
end

%scaling_factor = 1;
%power_factor = 1;
x = x.*scaling_factor;
x = x.^power_factor;
offset = min(min(x));
x = x-offset;


%A problem occurs if a potential lambda is negative and has absolute value larger than the smallest peak value.
% - this results in the glog algorithm trying to find sqrt of a negative number and producing complex results
% - therefore we should use fminbnd to search within a space constrained
% - around zero.

step_threshold = 1e-16;
options=optimset('TolX',step_threshold,'MaxFunEvals',1e3,'MaxIter',1e3);

%find the smallest value
small = min(x(:));
%lambda can be allowed to be negative only to the value of the smallest peak value squared. 
lower_limit = -small^2;
upper_limit = max(max(var(x'),(max(var(x'))/min(var(x'))))); %in some cases, the lambda value can be very large. Hopefully it is constrained by either of these options. 
lambda = fminbnd(@SSE, lower_limit,upper_limit,options, y0, x); 

lambda_std =5.0278*10^-09;
%TEST lambda for errors!!
error_flag = 0;
if abs(upper_limit-lambda)<=1e-5 %it seems that fminbind stops short of upper limit... 
    sprintf('glogopt error! lambda tending to infinity! Using standard lambda value %e\n', lambda_std);
    error_flag = 1;
elseif abs(lower_limit-lambda)<=1e-5
    sprintf('glogopt error! lambda tending to -infinity! Using standard lambda value %e\n', lambda_std);
    error_flag = 1;
end

% apply to all spectra for output

if isobject(dso) %if the QC subset was used, we will want to glog all samples
    x = dso.data;
    x = x';
    N = size(x, 2);
end

%replace optimised lambda with standard if errors have been recorded.
if error_flag
    lambda = lambda_std;
    sprintf('Standard lambda used - not optimal!!')
    %also scale the data - for whatever reason... ? (this was in original
    %normalise scripts etc) see Mark.
    
    %linear scaling
    %sum of intensities of all samples, average of that value
    %divide every sample by that average
    scaling_factor=sum(x);
    scaling_factor=mean(scaling_factor);
    scaling_factor = 1/scaling_factor;
    
    
end

x = x.*scaling_factor; %need to make all spectra scaled, not just QCs
x = x.^power_factor; %apply the power scaling if it has been adjusted...
x = x-min(x(:)); %make x sit on zero line.


%write these parameters to outfile
sprintf('Scaling Factor, %f\n',scaling_factor);
sprintf('Power Factor, %f\n',power_factor);
sprintf('Offset, %f\n',offset);
sprintf('Lambda, %f\n',lambda);
sprintf('Optimisation error, %i\n',error_flag);



for i=1:N
    z(:,i)=glog_dims(x(:,i),y0,lambda); %apply glog transform to scaled data using lambda
end

%% CREATE DATASET FOR OUTPUT
% Include all samples
include_out = 1:n;
dso.include{1} = include_out;

dso_out=dso;

if isobject(dso_out)
    dso_out.name=[dso.name,'_glogopt'];
    dso_out.data=z';
else
    dso_out = z';
end

%% SAVE OUTPUT DATASET
autoexport(dso_out, outfile, 'xml');

%Close the message output file
%fclose(fid);
 

%% RETURN PATH TO ORIGINAL STATE
if ~isdeployed
	path(original_path);
	rehash pathreset;
	rehash toolbox;
end

return
end

