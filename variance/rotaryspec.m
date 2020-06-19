%% function [myfreq,AUTOdata,ROTOdata,PHASE,COHER,slevel]=rotaryspec(specStruct)
%% Maarten Buijsman, IGPP, UCLA, 2008-02-13
%% computes auto, rotary, phase and coherence using 
%% uses input specStruct from functions cmgspecavg from CMGtool and routines out of cmgplotspectra  
%% output need to be used for plotting

function [myfreq,AUTOdata,ROTOdata,PHASE,COHER,slevel]=rotaryspec(specStruct);

myspec=specStruct.spec;
myfreq=specStruct.freq;
npieces=specStruct.npieces;
p=specStruct.conf;
nfft=specStruct.nfft;
nyqf = myfreq(end);

%% auto spectra
ydata=real(myspec(:,[1 2 3 4]));
[sm,sn]=size(ydata);
zz=min(real(min(ydata(:,1:sn-2)))); zz=log10(zz); zz=fix(zz); zz=10^zz;    
ydata=[ydata(:,1:sn-2) zz.*ydata(:,sn-1) zz.*ydata(:,sn)];
AUTOdata=ydata*nfft/(2*nyqf); %Convert Power spectra to Power spectral density     
AUTOdata = AUTOdata(:,1:2);

%% ROTARY spectra
ydata=real(myspec(:,[13 14 3 4]));
[sm,sn]=size(ydata);
zz=min(real(min(ydata(:,1:sn-2)))); zz=log10(zz); zz=fix(zz); zz=10^zz;    
ydata=[ydata(:,1:sn-2) zz.*ydata(:,sn-1) zz.*ydata(:,sn)];
ROTOdata=ydata*nfft/(2*nyqf); %Convert Power spectra to Power spectral density     
ROTOdata = ROTOdata(:,1:2);

%% phase
ydata=real(myspec(:,[11 11 11]));
[sm,sn]=size(ydata);
ydata=[ydata(:,1:sn-2) nan*ones(sm,1) nan*ones(sm,1)]; %% first two columns
nfftd=2; nyqfd=1;
PHASE=ydata*nfftd/(2*nyqfd); %Convert Power spectra to Power spectral density 
PHASE = PHASE(:,1);

%% coherence
ydata=real(myspec(:,[9 9 9]));
slevel=real(myspec(1,10));
[sm,sn]=size(ydata);
ydata=[ydata(:,1:sn-2) nan*ones(sm,1) nan*ones(sm,1)]; %% first two columns
nfftd=2; nyqfd=1;
COHER=ydata*nfftd/(2*nyqfd); %Convert Power spectra to Power spectral density 
slevel = slevel*ones(2,1);
COHER = COHER(:,1);