%% function R2 = varexp(x,y)
%  MCB, USM, 2021-3-21
%  computes how much of the variance in y (the fit) is explained by x
%  x should already have the mean removed
%  x and y are vectors and assumed equidistant in time (or space)
%  source: https://en.wikipedia.org/wiki/Coefficient_of_determination

function R2 = varexp(x,y);

% remove mean
y = y - mean(y);
x = x - mean(x);

% explained variance
R2 = ( 1 - mean((y-x).^2) / mean(y.^2) )*100;

