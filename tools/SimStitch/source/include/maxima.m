function [vals]=maxima(input_vec, number)
% maxima - return the locations of a given number of maxima 
%
% [output] = maxima(input, number of maxima)
%
% inputs : input_vec = a vector containing the input data
%          number = (obviously) the number of maxima to find
% 
% output : vals = vector containing the locations of the maxima found
%
% GPP, Uni of B'ham, 2001

%Do a quick check on the input incase the user's a fool
if ((size(input_vec, 1) ~= 1) & (size(input_vec, 2) ~= 1) | (ndims(input_vec) > 2))
    disp('Stop being a muppet. The input cannot be a matrix');
    return;
end

%Turn the input into a row-vector
if (size(input_vec, 1) > 1)
    input_vec = input_vec.';
end

%Firstly, find all the (local) maxima in the input
deriv = input_vec(2: size(input_vec,2)) - input_vec(1:size(input_vec,2)-1);
deriv_sq = deriv(2: size(deriv,2)) - deriv(1:size(deriv,2)-1);

mult = deriv(2: size(deriv,2)) .* deriv(1:size(deriv,2)-1);
turning_points = find(mult <= 0);

maxima = turning_points(deriv_sq(turning_points) < 0);

%Now, check if we have enough maxima - if not, just output those that we do have:
if (size(maxima,2) <= number)
    vals = maxima;
    return;
end

%Next, apply a value row above the maxima & sort the columns (ascending order)
maxima = sortrows([input_vec(maxima+1);maxima].').';

%Finally, output the locations of the biggest 'number' of rows
vals = fliplr(maxima(2,size(maxima,2)-number + 1:size(maxima,2)));
end

