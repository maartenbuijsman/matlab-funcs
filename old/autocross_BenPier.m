%% function [Gxx,fk,Nb] = autocross_BenPier(t,y1,y2,Kb,rem_ave,win);
%% Maarten Buijsman, NIOZ, 29-05-06
%% calculates autocorrelation and cross-correlation spectra according 
%% to Bendat and Piersol (1986), page 409
%% INPUT: vectors t (time) and y1 and y2 (data), Kb number of blocks, 
%% rem_ave to remove mean (not=0, yes=1) and win to indicate Hanning window
%% (1=yes) or no window smoothing (0=no)
%% OUTPUT: smoothed real-valued auto (y1 equals y2) or complex-valued cross-spectra 
%% (y1 is not y2) Gxx as a function of frequency fk, 
%% number of blocks Nb really used (2*Kb-1), 
%% and phase phi in radians (phi = 0 in case of autospectra), in time it can 
%% be expressed as DT=phi/(2*pi*fk); positive phi means that y2 leads y1; 
%% this gives positive DT and when added to y2 the data series must match

function [Gxx,fk,phi,Nb] = autocross_BenPier(t,y1,y2,Kb,rem_ave,win);

%% divide in Kb blocks and calculate spectra of y1 and y2
[X1,X2,Xh1,Xh2,fk,y1,y2,t,N,Ns,dt] = Xspectra_BenPier(t,y1,y2,Kb,rem_ave);

%% number of blocks really used
Nb = size(Xh1,1);

%% select yes/no window averaging
if win == 1; %% Hanning window
    Xh1 = Xh1; Xh2 = Xh2;    
else %% no window
    Xh1 = X1;  Xh2 = X2;    
end
    
%% 5) determine auto or cross-spectra and phase
Gxx = 2/(Nb*Ns*dt)*sum(conj(Xh1).*Xh2,1);
phi = atan2(-imag(Gxx),real(Gxx));

%% test
% Gii = 2/(Ns*dt)*(conj(Xh1).*Xh2);
% figure; loglog(fk,Gii); hold; loglog(fk,Gxx,'k.-');
% title(['win = ',num2str(win),'; Nb = ',num2str(Nb),'; rem ave = ',num2str(rem_ave)]);