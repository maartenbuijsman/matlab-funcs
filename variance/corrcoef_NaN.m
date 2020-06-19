%% function [r] = corrcoef_NaN(x,y)
%% Maarten Buijsman, NIOZ, 21-03-06
%% input: 2 simple vectors
%% removes NaNs and performs correlation analysis
%% returns r and propability p

function [r,p] = corrcoef_NaN(x,y);

%% make shure x and y are horizontal
if size(x,1)>1; x=x'; end
if size(y,1)>1; y=y'; end
    
%% remove NaNs
Inan1 = find(isnan(x));
Inan2 = find(isnan(y));
x([Inan1 Inan2]) = [];
y([Inan1 Inan2]) = [];

%% correlation
[RR,PP] = corrcoef(x,y); 
p = PP(2,1);
r = RR(2,1); 


