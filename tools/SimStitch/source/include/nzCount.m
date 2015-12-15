function nzCount = nzCount(x,dim);
% returns count of non-zero elements in x in given dimension

x(x~=0) = 1;
nzCount = sum(x,dim);
end

