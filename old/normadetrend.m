%% function [tn,xn] = normadetrend(t,x,step);
%% Maarten Buijsman, NIOZ, 04-08-2005
%% input time (t), variable (x), and time step (step)
%% step is used for regridding 
%% removes NaNs, grid, remove linear trend, and normalize with standar deviation
%% E&T, pg 329

function [tn,xn] = normadetrend(t,x,step);

%% remove NaNs
Inan = find(isnan(t));
t(Inan) = []; x(Inan) = [];
Inan = find(isnan(x));
t(Inan) = []; x(Inan) = [];

%% interp from beginning to end
% tn = [ceil(t(1)):step:floor(t(end))];
% xg = interp1(t,x,tn);

tn = t; xg=x; %% ungridded
[x_s,y_fit,r,p,cf] = line_fit(tn,xg,1);
xn = xg-(cf(1)*tn+cf(2));
xn = xn/std(xn);

