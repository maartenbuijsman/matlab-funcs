%% function [cflo,cfhi] = butter_freq_lowhigh(time,cutoff,N);
%% Maarten Buijsman, UCLA, 2010-03-05
%% works with bandpass2.m and bandpass3.m based on butter_freq.m
%% (gives identical results as butter_freq)
%% computes coefficients for low and high-passed solutions
%% INPUT:  vector time [days], cutoff that indicates period of the cutoff frequency, 
%% and N the order (1..9) (5 is OK)
%% OUTPUT: bf,af coefficients for low and high-passed filtering in structures cflo and cfhi

%% can also use butter(N,Wn,'stop') , which is the inverse of bandpass, 
%% however, it seems the frequency limits are inside the freq bands => removing more high freq stuff

function [cflo,cfhi] = butter_freq_lowhigh(time,cutoff,N);

cutoffreq    = 1./cutoff;               %% cycles per day
DT           = mean(diff(time)); 
if std(diff(time))>0.00001; disp('NOT EQUIDISTANT'); end
sampfreq     = 1/DT;                    
halfsampfreq = sampfreq*1/2;            %% Nyquist frequency
Wn           = cutoffreq/halfsampfreq;  %% normalize by Nyquist fr, see examples in Matlab
[cflo.bf,cflo.af] = butter(N,Wn,'low');
[cfhi.bf,cfhi.af] = butter(N,Wn,'high');
