function [d1,d2,mz2] = FillOut(d1,mz1,d2,mz2)
%make both inputs with same frequency points, extra data points are set to
%zero

%merge and sort
mz1 = [mz1 mz2];
[mz1 idx]=sortrows(mz1.');
mz1 = mz1.';

%sort data same way
d2 = [zeros(1,length(d1)) d2];
d1(length(mz1)) = 0;

d1 = d1(idx);
d2 = d2(idx);

%remove any duplicates
idx = find((mz1(1:(end-1))-mz1(2:end))==0);
if ~isempty(idx)
    mz1(idx) = [];
    d1(idx+1) = [];
    d2(idx) = [];
end

mz2 = mz1;
end

