%% function xf = butter_filt(x,bf,af);
%% Maarten Buijsman, UCLA, 2009-12-08
%% lowpass filtering using buttter, uses input from butter_freq.m
%% INPUT:  vector(!!!) x and coefficients [bf,af] from butter_freq
%% OUTPUT: xf

function xf = butter_filt(x,bf,af);

xf            = filtfilt(bf,af,x);
