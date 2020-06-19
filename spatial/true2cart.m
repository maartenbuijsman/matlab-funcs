%% [out] = true2cart(in)
%% Maarten Buijsman, NIOZ, 02-01-06
%% convert degrees positive clockwise relative to true north
%% to degrees cartesian (positive counterclockwise relative to x-axis)
%% in 0-180 and 180-0

function [out] = true2cart(in);

%%in = [-270:89:270];
%%in = [0:89:450];

%% map to cartesian
out = 90 - in;

%% get between 0 and 360
Isel = find(out>180);
out(Isel) = out(Isel)-360;

Isel = find(out<-180);
out(Isel) = out(Isel)+360;




