%%---------------------------------------------------------------------------------------------------
%%-- Function: rem_NaN
%%-- Author: Maarten Buijsman
%%-- Place: GFDL
%%-- Date: 6-1-2011
%%-- Description: 
%%-- Function removes NaNs
%%---------------------------------------------------------------------------------------------------
function Xout = rem_NaN(X);
%Inotnan = find(~isnan(X));
Inotnan = find(~isnan(X) & ~isinf(abs(X))); %also remove inf values
Xout = X(Inotnan);

