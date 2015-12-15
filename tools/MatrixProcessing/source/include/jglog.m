function [zj, gmn]=jglog(y,y0,lambda) % Rescale variables using Jacobian: w = J*z
% Note slight difference in format to Durbin paper - makes eqn computational
% (has extra multiplicative term only; moving minimum up)

z=glog(y,y0,lambda);
gmn=exp(mean(log(sqrt((y-y0).^2+lambda))));
zj=z.*gmn;
end
