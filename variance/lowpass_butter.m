%% function [xf] = lowpass_butter(time,x,cutoff,N);
%% Maarten Buijsman, UCLA, 2009-04-21
%% lowpass filtering using buttter
%% INPUT:  vectors time [days] and x, and cutoff, which indicates period of the cutoff frequency 
%% N the order (1..9) (5 is OK)
%% OUTPUT: lowpass filtered vector xf

function [xf] = lowpass_butter(time,x,cutoff,N);

cutoffreq    = 1./cutoff;               %% cycles per day
DT           = mean(diff(time)); 
sampfreq     = 1/DT;                    
halfsampfreq = sampfreq*1/2;            %% Nyquist frequency
Wn           = cutoffreq/halfsampfreq;  %% normalize by Nyquist fr, see examples in Matlab
%Wn = interp1([0 halfsampfreq],[0 1],1/2*cutoffreq); %% INCORRECT!!
[bf,af]      = butter(N,Wn,'low');
xf           = filtfilt(bf,af,x);

