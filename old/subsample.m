%% function varout = subsample(tin,varin,tbase)
%% Maarten Buijsman, IGPP, UCLA, 2010-05-12
%% subsamples by computing the average value over
%% a time step of the new time base tbase
%% all time series need to be equidistant
%% NOT SURE WHAT HAPPENS WHEN tbase is larger then tin

function varout = subsample(tin,varin,tbase);

% test
% tin = ti;
% varin = ug(i,:);
% tbase = tt;

dt        = tbase(2)-tbase(1);
varinsum  = cumtrapz(tin,varin);

% %% when gradient is used
% varinsumi = interp1(tin,varinsum,tbase);  
% varout    = gradient(varinsumi,dt);       

%% alternate method is better
tbase(end+1) = tbase(end)+dt;                     %% extend time base to acommodate shift in time
varinsumi    = interp1(tin,varinsum,tbase-dt/2,'linear','extrap');  %% goes wrong when tbase is larger; make safeguard
%varinsumi    = interp1(tin,varinsum,tbase-dt/2);  %% shift in time
varout       = diff(varinsumi)/dt;                %% to get output on tbase

return



% figure
% plot(tin,varin)
% hold
% plot(tbase(1:end-1),varout,'r-')