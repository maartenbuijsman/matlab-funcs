%% function [uh,us,ud,ul] = bandpass_notail.m(utot,bf0,af0,bfs,afs,bfd,afd);
%% Maarten Buijsman, GFDL, 2011-05-06
%% INPUT equidistant time series band utot 
%% pass filter using coefficients from butter_freq.m
%% bf0,af0 (upper bound high-pass),
%% bfs,afs (upper bound semidiurnal),
%% bfd,afd (upper bound diurnal)
%% OUTPUT high-passed uh, s.d. passed us, d. passed ud, low-passed ul
%% 2010-01-17: added mirrored head and tail to reduce ringing (nothing is ideal....)

function [uh,us,ud,ul] = bandpass_notail(ut,bf0,af0,bfs,afs,bfd,afd);

%ut=utot(1,:)

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

%% filter --------------------------------------
xf0 = filtfilt(bf0,af0,ut); uh = ut - xf0;           % high-passed frequencies
xfs = filtfilt(bfs,afs,ut); us = ut - xfs - uh;      % semidiurnal frequencies
ul  = filtfilt(bfd,afd,ut); ud = ut - ul  - us - uh; % diurnal frequencies

% %% remove head and tail --------------------------------------
% uh=uh(Idum);
% us=us(Idum);
% ud=ud(Idum);
% ul=ul(Idum);

% figure
% plot(ud,'r-')
% hold
% plot(ut,'k-')

%test (OK)
%ut = uh + us + ud + ul; ut(1:10)
%utot(1:10)