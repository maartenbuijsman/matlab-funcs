%% function [bf,af] = butter_freq_band(time,cutoffhi,cutofflo,N);
%% Maarten Buijsman, UCLA, 2010-03-05
%% works with bandpass3.m based on butter_freq.m
%% computes coefficients for low and high-passed solutions
%% INPUT:  vector time [days], cutoffhi and cutofflo
%% that indicate PERIOD [days] of the low and high cutoff frequency, 
%% and N the order (1..9) (5 is OK); with butter use higher order to get steeper drop-off
%% OUTPUT: bf,af coefficients for band-passed filtering

function [bf,af] = butter_freq_band(time,cutofflo,cutoffhi,N);

DT           = mean(diff(time)); 
if std(diff(time))>0.00001; disp('NOT EQUIDISTANT'); end
sampfreq     = 1/DT;                    
halfsampfreq = sampfreq*1/2;            %% Nyquist frequency

%% gives low freq
cutoffreqhi    = 1./cutoffhi;               %% cycles per day
Wnhi           = cutoffreqhi/halfsampfreq;  %% normalize by Nyquist fr, see examples in Matlab

%% gives high freq
cutoffreqlo    = 1./cutofflo;               %% cycles per day
Wnlo           = cutoffreqlo/halfsampfreq;  %% normalize by Nyquist fr, see examples in Matlab

[bf,af] = butter(N,[Wnhi Wnlo],'bandpass');

