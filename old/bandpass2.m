%% function [uh,us,ud,ul] = bandpass2(utot,bf0,af0,bfs,afs,bfd,afd);
%% Maarten Buijsman, UCLA, 2010-03-05
%% takes input from butter_freq_lowhigh.m
%% (gives identical results as bandpass)
%% INPUT equidistant time series band utot 
%% bf0,af0 (upper bound high-pass),
%% bfs,afs (upper bound semidiurnal),
%% bfd,afd (upper bound diurnal)
%% OUTPUT high-passed uh, s.d. passed us, d. passed ud, low-passed ul
%% 2010-01-17: added mirrored head and tail to reduce ringing (nothing is ideal....)

function [uh,us,ud,ul] = bandpass2(ut,cflow0,cfhi0,cflows,cfhis,cflowd,cfhid);

%% add head and tail 10%  --------------------------------------
le = length(ut);
Iextra = floor(le*10/100);
uhead = ut(Iextra:-1:1);
utail = ut(end:-1:end-Iextra+1);

Idum = Iextra+1:Iextra+le;
udum(1:Iextra)=uhead;
udum(Idum)=ut;
udum(Iextra+le+1:2*Iextra+le)=utail;
ut = udum;
%size(udum),size(uhead),size(utail)

%% filter --------------------------------------
uh = filtfilt(cfhi0.bf,cfhi0.af,ut); %% high freq stuff

u1 = filtfilt(cflow0.bf,cflow0.af,ut); us = filtfilt(cfhis.bf,cfhis.af,u1);  %% semiudiurnal

u1 = filtfilt(cflows.bf,cflows.af,ut); ud = filtfilt(cfhid.bf,cfhid.af,u1);  %% diurnal

ul = filtfilt(cflowd.bf,cflowd.af,ut); %% subtidal stuff


%% remove head and tail --------------------------------------
uh=uh(Idum);
us=us(Idum);
ud=ud(Idum);
ul=ul(Idum);

