%% function [uh,us,ud,ul] = bandpass3(ut,cfhi0,bfs,afs,bfd,afd,cflod);
%% Maarten Buijsman, UCLA, 2010-03-05
%% takes input from butter_freq_lowhigh.m and butter_freq_band 
%% INPUT equidistant time series band utot 
%% OUTPUT high-passed uh, s.d. passed us, d. passed ud, low-passed ul
%% 2010-01-17: added mirrored head and tail to reduce ringing (nothing is ideal....)

function [uh,us,ud,ul] = bandpass3(ut,cfhi0,bfs,afs,bfd,afd,cflod);

% %% add head and tail 10%  --------------------------------------
% le = length(ut);
% Iextra = floor(le*10/100);
% uhead = ut(Iextra:-1:1);
% utail = ut(end:-1:end-Iextra+1);
% 
% Idum = Iextra+1:Iextra+le;
% udum(1:Iextra)=uhead;
% udum(Idum)=ut;
% udum(Iextra+le+1:2*Iextra+le)=utail;
% ut = udum;
% %size(udum),size(uhead),size(utail)

%% filter (gives identical results as bandpass) --------------------------------------
uh = filtfilt(cfhi0.bf,cfhi0.af,ut); %% high freq stuff

us = filtfilt(bfs,afs,ut);  %% semiudiurnal

ud = filtfilt(bfd,afd,ut);  %% diurnal

ul = filtfilt(cflod.bf,cflod.af,ut); %% subtidal stuff


% %% remove head and tail --------------------------------------
% uh=uh(Idum);
% us=us(Idum);
% ud=ud(Idum);
% ul=ul(Idum);

