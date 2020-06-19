%% function [bf,af] = butter_freq(time,cutoff,N);
%% Maarten Buijsman, UCLA, 2009-12-08
%% determines frequency information, coefficients to be used in butter_filt.m and bandpass.m
%% INPUT:  vector time [days], cutoff that indicates period of the cutoff frequency, 
%% and N the order (1..9) (5 is OK)
%% OUTPUT: bf,af coefficients to be used in butter_filt.m

function [bf,af] = butter_freq(time,cutoff,N);

cutoffreq    = 1./cutoff;               %% cycles per day
%DT           = mean(diff(time)); 
DT           = time(2)-time(1);
if std(diff(time))>0.00001; disp('NOT EQUIDISTANT'); end
sampfreq     = 1/DT;                    
halfsampfreq = sampfreq*1/2;            %% Nyquist frequency
Wn           = cutoffreq/halfsampfreq;  %% normalize by Nyquist fr, see examples in Matlab
[bf,af]      = butter(N,Wn,'low');
%xf           = filtfilt(bf,af,x);
