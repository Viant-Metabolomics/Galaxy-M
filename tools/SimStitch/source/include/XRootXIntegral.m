function [myInt] = XRootXIntegral(z,k,a,b,c,x,q)
% determines the integral of the function R^z.sqrt(R)
% at the limit x 
% where R = a*x^2 + b*x + c
% and q = 4*a*c - b^2
% and k = 4*a/q
% 
% NOTE: code optimised for x specified at y=0 values: will not work for
% other values of x!!

if a>0
    inta = 1/sqrt(a)*asinh((2*a*x+b)/sqrt(q));
else
    inta = -1/sqrt(-a)*asin((2*a*x+b)/sqrt(-q));
end
R = (a*x^2 + b*x + c);
%     suma = sum(factorial([0:z]) .* factorial([1:z+1]) .* (4*k*R).^([0:z]) ./ factorial(2*[0:z]+2));
% note that since we are looking for the zero-crossing, R will always
% be tending ~0, and so suma ~ 0.5
suma = 0.5;
myInt = factorial(2*z + 2) / ((factorial(z+1))^2 * (4*k)^(z+1)) * ( k*(2*a*x + b)*sqrt(R) / a * suma + inta );
end

