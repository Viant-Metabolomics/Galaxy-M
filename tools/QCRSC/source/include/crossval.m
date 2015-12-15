function loss = crossval(funorStr,varargin)
%CROSSVAL cross-validation
%   VALS = CROSSVAL(FUN,X) performs 10-fold cross-validation for the
%   function FUN, applied to the data in X. FUN is a function handle that
%   is called 10 times with two inputs, the training set of X and the test
%   set of X, as follows:
%
%      TESTVAL = FUN(XTRAIN,XTEST)
%
%   Each time it is called, FUN should use XTRAIN to fit a model, then
%   return some criterion TESTVAL computed on XTEST using that fitted
%   model. X can be a column vector or a matrix. Rows of X correspond to
%   observations; columns correspond to variable or features. VALS can be a
%   column vector or a matrix. Each row of VALS contains the result of
%   applying FUN to one test set. If FUN returns a non-scalar value,
%   CROSSVAL reshapes it to a row vector using linear indexing order and
%   stores it in a row of the matrix VALS.
%
%   VALS = CROSSVAL(FUN,X,Y,...) is used when there are multiple variables
%   X, Y, ....  All variables (column vectors, matrices or arrays) must
%   have the same number of rows. FUN is a function handle that is called
%   with the training subsets of X, Y, ..., followed by the test subsets of
%   X, Y, ..., as shown below:
%
%      TESTVALS = FUN(XTRAIN,YTRAIN,...,XTEST,YTEST,...)
%
%   MSE = CROSSVAL('mse',X,Y,'Predfun', PREDFUN) returns MSE, a scalar
%   containing a 10-fold cross-validation estimate of mean-squared error
%   for the function PREDFUN. X is the predictor which can be a column
%   vector, a matrix or an array. Y is a numeric column vector containing
%   the response values. X and Y must have the same number of rows. PREDFUN
%   is a function handle that is called with the training subset of X,
%   followed by the training subset of Y, and then the test subset of X as
%   follows:
%
%      YFIT = PREDFUN(XTRAIN,YTRAIN,XTEST)
%
%   Each time it is called, PREDFUN should use XTRAIN and YTRAIN to fit
%   a regression model and then return fitted values in a column vector
%   YFIT. Each row of YFIT contains the predicted values for the
%   corresponding row of XTEST. CROSSVAL computes the squared errors
%   between YFIT and the corresponding response test set, and returns the
%   overall mean across all test sets.
%
%   MCR = CROSSVAL('mcr',X,Y,'Predfun', PREDFUN) returns MCR, a scalar
%   containing a 10-fold cross-validation estimate of misclassification
%   rate (the proportion of misclassified samples) for the function PREDFUN
%   with the matrix X as predictor values and vector Y as class labels.
%   PREDFUN should use XTRAIN and YTRAIN to fit a classification model and
%   return YFIT as the predicted class labels for XTEST. CROSSVAL then
%   computes the number of misclassifications between YFIT and the
%   corresponding response test set, and returns the overall
%   misclassification rate across all test sets.
%
%   VAL = CROSSVAL(CRITERION,X1,X2,...,Y,'Predfun', PREDFUN) ), where
%   CRITERION is 'mse' or 'mcr', returns a cross-validation estimate of
%   mean-squared error (for a regression model) or misclassification rate
%   (for a classification model) with X1, X2, ... as predictor values and
%   Y as response values or class labels. All variables (X1, X2, ... Y)
%   must have the same number of rows. PREDFUN is a function handle that is
%   called with the training subsets of X1, X2, ..., followed by the
%   training subset of Y, and then the test subsets of X1, X2, ..., as
%   shown below:
%
%      YFIT = PREDFUN(X1TRAIN,X2TRAIN,...,YTRAIN,X1TEST,X2TEST,...)
%
%   YFIT should be a column vector containing the fitted values.
%
%   VALS = crossval(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs chosen from the following:
%
%   'Kfold'        The number of folds K for K-fold cross-validation.
%   'Holdout'      The ratio or the number of observations P for holdout.
%                  P must be a scalar. When 0<P<1, it randomly selects
%                  approximately P*N observations for the test set. When P is
%                  an integer, it randomly selects P observations for the
%                  test set.
%   'Leaveout'     The value must be 1. Leave-one-out cross-validation is
%                  used.
%   'Partition'    A CVPARTITION object C.
%   'Stratify'     A column vector GROUP indicating the group information
%                  for stratification. Both training and test sets have
%                  roughly the same class proportions as in GROUP.
%                  CROSSVAL treats NaNs or empty strings in GROUP as
%                  missing values and the corresponding rows of X,Y,...
%                  are ignored. A stratified partition is preferred for
%                  evaluating classification algorithms.
%   'Mcreps'       A positive integer indicating the number of Monte-Carlo
%                  repetitions for validation. If the first input of
%                  CROSSVAL is 'mse' or 'mcr', CROSSVAL returns the mean of
%                  mean-squared error or misclassification rate across all
%                  of the Monte-Carlo repetitions. Otherwise, CROSSVAL
%                  concatenates the values of VALS from all the Monte-Carlo
%                  repetitions along the first dimension.
%   'Options'      A structure that contains options specifying whether to
%                  conduct multiple function evaluations in parallel, and
%                  options specifying how to use random numbers when computing
%                  cross validation partitions. This argument can be created
%                  by a call to STATSET. CROSSVAL uses the following fields:
%                      'UseParallel'
%                      'UseSubstreams'
%                      'Streams'
%                  For information on these fields see PARALLELSTATS.
%                  NOTE: If supplied, 'Streams' must be of length one.
%
%   Only one of 'Kfold', 'Holdout','Leaveout', or 'Partition' parameters
%   can be specified.  10-fold cross-validation is the default if none of
%   the 'Kfold', 'Holdout', 'Leaveout' or 'Partition' parameters is
%   provided. The 'Partition' parameter cannot be specified with
%   'Stratify'. If the values of both 'Partition' and 'Mcreps' are
%   provided, the first Monte-Carlo repetition will use the partition
%   information contained in the given CVPARTITION object, the REPARTITION
%   method will be called on this CVPARTITION object to generate a new
%   partition for each of the remaining Monte-Carlo repetitions.
%
%   Examples:
%      % Compute mean-square error for regression using 10-fold
%      % cross-validation.
%      load('fisheriris');
%      y = meas(:,1);
%      x = [ones(size(y,1),1) meas(:,2:4)];
%      regf = @(xtrain, ytrain,xtest)(xtest * regress(ytrain,xtrain));
%      cvMse = crossval('mse',x,y,'predfun',regf)
%
%      % Compute misclassification rate using stratified 10-fold
%      % cross-validation
%      load('fisheriris');
%      y = species;
%      % A stratified partition is preferred to evaluate classification
%      % algorithms.
%      cp = cvpartition(y,'k',10);
%      classf = @(xtrain, ytrain,xtest)(classify(xtest,xtrain,ytrain));
%      cvMCR = crossval('mcr',meas,y,'predfun', classf,'partition',cp)
%
%      % Return the confusion matrix using stratified 10-fold
%      % cross-validation.
%      load('fisheriris');
%      % Each of the ten confusion matrices needs to sort the group labels
%      % according to the same order
%      yorder = unique(species);
%      % A stratified partition is preferred to evaluate classification
%      % algorithms.
%      cp = cvpartition(species,'k',10);
%      f = @(xtr,ytr,xte,yte) confusionmat(yte,classify(xte,xtr,ytr),...
%          'order',yorder);
%      cfMat = crossval(f,meas,species,'partition',cp);
%      % cfMat is the summation of 10 confusion matrices from 10 test sets.
%      cfMat = reshape(sum(cfMat),3,3)
%
%   See also CVPARTITION, CVPARTITION/REPARTITION, STATSET, PARALLELSTATS.

%   References:
%   [1] Hastie, T. Tibshirani, R, and Friedman, J. (2001) The Elements of
%       Statistical Learning, Springer, pp. 214-216.

% Copyright 2007-2012 The MathWorks, Inc.


if nargin < 2
    error(message('stats:crossval:TooFewInputs'));
end

firstInputType = 'func';
ismse='';
if ischar(funorStr)
    funorStr = lower(funorStr);
    if strcmp(funorStr,'mse')
        ismse = true;
        firstInputType = 'lossMeasure';
    elseif strcmp(funorStr,'mcr')
        ismse = false;
        firstInputType = 'lossMeasure';
    else
        error(message('stats:crossval:BadFun'));
    end
elseif ~isa(funorStr,'function_handle')
    error(message('stats:crossval:BadFun'));
end

n = size(varargin{1},1);
if n < 2
    error(message('stats:crossval:TooFewDataRows'));
end

nData = length(varargin);
for i = 2:length(varargin);
    if size(varargin{i},1) ~= n
        if ~(ischar(varargin{i}) && size(varargin{i},1) ==1 )
            error(message('stats:crossval:MismatchedDataRows'));
        else
            nData = i-1;
            break;
        end
    end
end

if strncmp(firstInputType, 'lossMeasure',3) && nData < 2
    error(message('stats:crossval:NotEnoughVars'));
end

data = varargin(1:nData);
varargin(1:nData)= [];

pnames = { 'kfold'  'holdout'  'leaveout' 'mcreps' 'stratify' 'partition' 'predfun' ...
    'options' };
dflts =  { []        []          []        1          []         []       []        ...
    statset('crossval') };
[nfolds,holdout, leaveout, mcreps, stratify, cvp, predfun, parallelOption] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

if strncmp(firstInputType,'lossMeasure',3)
    
    if isempty(predfun)
        error(message('stats:crossval:MissingPredfun'));
    end
       
    if ~ismse %mcr
        if ischar(data{end})  && ismatrix(data{end})
            data{end}=nominal(data{end});
        elseif (iscellstr(data{end}) || isnumeric(data{end}) || ...
                islogical(data{end}) || isa(data{end},'categorical') ) && ...
                isvector(data{end})
            data{end}=nominal(data{end});
        else
            error(message('stats:crossval:MCRBadGroup'));
        end
    else %mse
        if ~isfloat(data{end}) || ~isvector(data{end})
            error(message('stats:crossval:MSEBadY'));
        end     
    end
    
else %the first input is a function handle
    if ~isempty(predfun)
        warning(message('stats:crossval:UnneededPredfun'));
    end
end

if ~isempty(leaveout) && leaveout ~= 1
    error(message('stats:crossval:UnsupportedLeaveout'));
end

if ~( isnumeric(mcreps) && isscalar(mcreps) && mcreps == round(mcreps)...
        && mcreps >= 1)
    error(message('stats:crossval:BadMcreps'));
end

choices=[];

cvopts = sum( [~isempty(holdout), ~isempty(nfolds) , ~isempty(leaveout)]);
if cvopts > 1
    error(message('stats:crossval:InconsistentOpts'));
elseif ~isempty(cvp)
    if cvopts > 0 || ~isempty(stratify)
        error(message('stats:crossval:ConflictOpts'));
    elseif ~isa(cvp,'cvpartition')
        error(message('stats:crossval:Badcvp'));
    elseif cvp.N ~= n
        error(message('stats:crossval:MissizedCVPN'));
    end
    choices = 'cvpartition';
    if mcreps > 1 && (strcmp(cvp.Type,'resubstitution') ...
            || strcmp(cvp.Type,'leaveout'))
        mcreps = 1;
        warning(message('stats:crossval:InvalidMcreps',cvp.Type));
    end
else
    if ~isempty(stratify)
        stratify = grp2idx(stratify);
        if size( stratify,1) ~= n
            error(message('stats:crossval:MissizedStratify'));
        end
        
    end
    if cvopts == 0
        choices = 'kfold';
        cvarg = 10; %default cross-validation option
    elseif ~isempty(nfolds)
        choices = 'kfold';
        cvarg = nfolds;
    elseif ~isempty(holdout)
        choices ='holdout';
        cvarg = holdout;
    elseif ~isempty(leaveout)
        choices = 'leaveout';
        if mcreps > 1
            mcreps = 1;
            warning(message('stats:crossval:InvalidMcreps',choices));
       end
        cvarg = 1;
    end
end

% Process options parallel computation and random number streams.
[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions( ...
        parallelOption,false);
useSubstreams    = RNGscheme.useSubstreams;
streams          = RNGscheme.streams;
useDefaultStream = RNGscheme.useDefaultStream;

if useDefaultStream
    s = RandStream.getGlobalStream;
else
    s = streams{1};
end

substreamOffset = internal.stats.parallel.freshSubstream(s);
substreamOffset = substreamOffset -1 ;

% Process cvpartition(s)
if ~strcmp(choices,'cvpartition')
    if isempty(stratify)
        cvp = cvpartition(n, choices,cvarg, s);
    else
        cvp = cvpartition(stratify, choices,cvarg, s);
    end
    % We used up a Substream to make the first cvpartition.
    substreamOffset = substreamOffset+1;
end

nData = size(data,2);

numtests = cvp.NumTestSets;
numiter = numtests * mcreps;

if useParallel
    % Get all re-randomized cv partitions up front
    cvparray = cell(mcreps,1);
    cvparray{1} = cvp;
    for j=2:mcreps
        if useSubstreams
            set(s,'Substream',j+substreamOffset);
        end
        cvparray{j} = cvp.repartition(s);
    end
end

% Perform the function evaluations
if strcmp(firstInputType,'func') % the first input is a function handle
    
    % Initialize before looping to allocate the output array
    funResult = getFuncVal(1, nData, cvp, data, funorStr, []);
    szFunResult = size(funResult);
    loss(numiter,:) = funResult(:)';
    loss(1,:) = funResult(:)';
    
    if useParallel
        parfor j=2:numiter
            mcrep = 1 + floor((j-1)/numtests);
            i = 1 + mod(j-1,numtests);  % i-th fold within this mc rep
            cvpj = cvparray{mcrep};
            funResult = getFuncVal(i, nData, cvpj, data, funorStr, szFunResult);
            
            % Accumulate the results
            loss(j,:) = funResult(:)';
        end  % parfor mr=1:mcreps
    else
        for mr = 1:mcreps
            if (mr==1)
                loopstart=2;
            else
                loopstart=1;
                if useSubstreams
                    set(s,'Substream',mr+substreamOffset);
                end
                cvp = cvp.repartition(s);
            end
            
            offset = cvp.NumTestSets * (mr-1);
            
            for i = loopstart:cvp.NumTestSets
                funResult = getFuncVal(i, nData, cvp, data, funorStr, szFunResult);
                loss(offset+i,:) = funResult(:)';
            end
        end  % for mr=1:mcreps
    end
    
else % the first input is the loss measure
    
    loss = zeros(numtests,1);
          
    if useParallel
        
        parfor j=1:numiter
            mcrep = 1 + floor((j-1)/numtests);
            i = 1 + mod(j-1,numtests);  % i-th fold within this mc rep
            cvpj = cvparray{mcrep};
            
            [funResult,outarg] = getLossVal(i, nData, cvpj, data, predfun);
            
            % Accumulate the results
            if ismse
                temploss = sum((outarg-funResult).^2);
            else
                if ~isa(funResult, 'nominal')
                    funResult = nominal(funResult);
                end
                temploss = sum(outarg ~= funResult);
            end
            loss(j,:) = temploss;
        end
        
    else
        
        for mr = 1:mcreps
            if mr > 1
                if useSubstreams
                    set(s,'Substream',mr+substreamOffset);
                end
                cvp = cvp.repartition(s);
            end
            offset = cvp.NumTestSets * (mr-1);
            
            for i = 1:cvp.NumTestSets
                [funResult,outarg] = getLossVal(i, nData, cvp, data, predfun);
                
                % Accumulate the results
                if ismse
                    temploss = sum((outarg-funResult).^2);
                else
                    if ~isa(funResult, 'nominal')
                        funResult = nominal(funResult);
                    end
                    temploss = sum(outarg ~= funResult);
                end
                loss(offset+i,:) = temploss;
            end
        end
        
    end
    
    %loss = sum(loss)/ (mcreps * sum(cvp.TestSize));
    loss = median(loss);
    
end  % if/else firstInputType

%
% === Clean-up ====
%

% If the user passed a 'Streams' argument, then we may need to update
% or restore state of the stream(s) that were passed.
% If neither 'UseParallel' nor 'UseSubstreams' was selected, however,
% we do nothing here, because we want the effects we had on the stream state
% to persist outside the function call.

if useSubstreams
    
    % If we are here, the user passed in a stream, which we used.
    % On entry, that stream had Substream set to substreamOffset.
    % If no cvpartition was passed to us, we constructed one, and
    % incremented substreamOffset to reflect that we had "used up"
    % one Substream.  Afterwards, we used an additional Substream
    % for each of the remaining (mcreps-1) monte carlo repetitions.
    % On exit, we increment the Substream just beyond the last value
    % that we used. This will keep sequences produced by the stream
    % within this function segregated (ie, non-overlapping) from
    % streams produced outside the function, assuming that outside
    % the function the user does not reassign the Substream within
    % the range that we used. This action occurs on the client.
    % There is no need (or advantage) in going to the workers,
    % even if parallel.
    
    %s = assignStream(streams,useSubstreams);
    set(s,'Substream',substreamOffset+mcreps);
    
end

end  % of crossval()


function funResult = evalFun(fun,arg)
try
    funResult = fun(arg{:});
catch ME
    if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
            && ~isempty(strfind(ME.message, func2str(fun)))
        error(message('stats:crossval:FunNotFound', func2str( fun )));
    else
        error(message('stats:crossval:FunError', func2str( fun ), ME.message));
    end
end

end  % of evalFun()

function funResult = getFuncVal(i, nData, cvp, data, funorStr, szFunResult)
arg = cell(2*nData,1);
train = cvp.training(i);
test = cvp.test(i);
% Take subsets of the inputs
for k = 1:nData
    arg{k} = data{k}(train,:);
    arg{nData+k} = data{k}(test,:);
end

% Apply the function to the current subset
funResult = evalFun(funorStr,arg(:));

% Check that size is okay
if ~isempty(szFunResult) && ~isequal(size(funResult),szFunResult)
    error(message('stats:crossval:FunOutSizeMismatched', func2str( funorStr )));
end

end  % of funResult()

function [funResult,outarg] = getLossVal(i, nData, cvp, data, funorStr)
arg = cell(2*nData,1);
train = cvp.training(i);
test = cvp.test(i);
% Take subsets of the inputs
for k = 1:nData
    arg{k} = data{k}(train,:);
    arg{nData+k} = data{k}(test,:);
end

% Apply the function to the current subset
funResult = evalFun(funorStr,arg(1:end-1));
outarg = arg{end};

% Check that the size is okay
if ~isequal(size(funResult), size(outarg) )
    error(message('stats:crossval:badFunResultSize'));
end

end  % of getLossVal()
