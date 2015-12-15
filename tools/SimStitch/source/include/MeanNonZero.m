function y = MeanNonZero(x,dim)
%MEAN   Average or mean value.
% 
% 23/Oct/07 TGP     Mean of only non-zero elements
%   For vectors, MEAN(X) is the mean value of the elements in X. For
%   matrices, MEAN(X) is a row vector containing the mean value of
%   each column.  For N-D arrays, MEAN(X) is the mean value of the
%   elements along the first non-singleton dimension of X.
%
%   MEAN(X,DIM) takes the mean along the dimension DIM of X. 
%
%   Example: If X = [0 1 2
%                    3 4 5]
%
%   then mean(X,1) is [1.5 2.5 3.5] and mean(X,2) is [1
%                                                     4]
%
%   Class support for input X:
%      float: double, single
%
%   See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.
%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.17.4.3 $  $Date: 2005/05/31 16:30:46 $

if nargin==1, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end

  nzCnt = nzCount(x,dim);
  nzCnt(nzCnt==0) = eps;    %avoid div by 0 errors
  y = sum(x)./nzCnt;
else
  nzCnt = nzCount(x,dim);
  nzCnt(nzCnt==0) = eps;
  y = sum(x,dim)./nzCnt;
end

end

