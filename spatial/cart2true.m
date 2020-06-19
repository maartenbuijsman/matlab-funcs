%% [out] = cart2true(in)
%% Maarten Buijsman, NIOZ, 02-01-06
%% convert degrees cartesian (positive counterclockwise relative to x-axis)
%% to degrees positive clockwise relative to true north

function [out] = cart2true(in);

%%in = [-270:90:270];

%% map to true
out = 90 - in;

%% get between 0 and 360
Isel = find(out>360);
out(Isel) = out(Isel)-360;

Isel = find(out<0);
out(Isel) = out(Isel)+360;




