%% function [gam2,phi,fk,yam2_conf,EDOF,Nb] = coherence_BenPier(t,y1,y2,Kb,conf,rem_ave,win);
%% Maarten Buijsman, NIOZ, 25-05-06
%% calculates coherence according to Bendat and Piersol (1986), page 409
%% INPUT: vectors t (time) and y1 and y2 (data), Kb number of blocks, 
%% %confidence level conf (from E&T p488; but is not yet fully understood...) 
%% rem_ave to remove mean (not=0, yes=1) and 
%% win to indicate Hanning window (1=yes) or no window smoothing (0=no)
%% note that t=n*dt, n=0,1,2,...N-1, N is length data
%% OUTPUT: coherence squared gam2, phase in radians phi, 
%% frequencies fk=k/(N*dt) n=0,1,2,...N/2, fn = Nyquist frequency (this is the cutoff frequency),
%% %-conf confidence level yam2_conf, equivalent degrees of freedom EDOF
%% and number of blocks Nb really used (2*Kb-1)
%% NOTE: positive phi: y2 leads y1
%%       negative phi: y1 leads y2 (see test below)

function [gam2,phi,fk,yam2_conf,EDOF,Nb] = coherence_BenPier(t,y1,y2,Kb,conf,rem_ave,win);

%% divide in Kb blocks and calculate spectra of y1 and y2
[X1,X2,Xh1,Xh2,fk,y1,y2,t,N,Ns,dt] = Xspectra_BenPier(t,y1,y2,Kb,rem_ave);

%% select yes/no window averaging
if win == 1; %% Hanning window
    Xh1 = Xh1; Xh2 = Xh2;    
else %% no window
    Xh1 = X1;  Xh2 = X2;    
end

%% number of blocks
Nb = size(Xh1,1);

%% autospectral power density 
G11 = 2/(Nb*Ns*dt)*sum(conj(Xh1).*Xh1,1);
G22 = 2/(Nb*Ns*dt)*sum(conj(Xh2).*Xh2,1);
%figure; loglog(fk,G22); hold; plot([1 1]*fqz,[1e0 1e2],'r--')

%% cross-spectral power density and average over blocks
G12 = 2/(Nb*Ns*dt)*sum(conj(Xh1).*Xh2,1);

%% coherence
gam2 = abs(G12).^2./(G11.*G22);

%% phase B&P, p408 
%% note that -j => clockwise rotating circle (!!!)
%% positive phi means that y2 leads y1
%% negative phi means that y1 leads y2 (see test below)
phi  = real(log(G12./abs(G12))/(-j));
%phib = atan2(-imag(G12),real(G12));

%% get uncertainty
%% use Hanning EDOF = 8/3*N/(1/2*Ns)*Nb ?? without Nb
%EDOF = 8/3*N/(1/2*Ns);
EDOF = Nb;
conf = 95;
yam2_conf = 1-((100-conf)/100)^(1/(EDOF-1));

%% plot
% figure
% subplot(2,1,1)
% plot(fk,gam2); hold
% plot(fk,yam2_conf*ones(size(gam2)),'k--'); 
% ylabel('\gamma_{12}^2')
% 
% subplot(2,1,2)
% plot(fk,phi*180/pi)
% xlabel('fk [Hz]')
% ylabel('\eta [^o]')

%% test coherence phase shift
%% y1 leads y2 => negative phase
% t = [0:1:12*40]; y1 = cos(2*pi/12*t); y2 = 0.8*cos(2*pi/12*t-30*pi/180); figure; plot(t,y1,'r-',t,y2,'b-'); title('red leads blue')
% Kb = 5;conf =95;rem_ave=1;,win=1;
% [gam2,phi,fk,yam2_conf,EDOF,Nb] = coherence_BenPier(t,y1,y2,Kb,conf,rem_ave,win);
% figure; subplot(2,1,1); plot(fk,gam2); subplot(2,1,2); plot(fk,phi*180/pi);



