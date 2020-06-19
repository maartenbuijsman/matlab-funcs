%% function [Igap] = find_gap(time,time_gap);
%% Maarten Buijsman, NIOZ, 28-09-05
%% function find gaps in time larger than time_gap
function [Igap] = find_gap(time,time_gap);

% time_gap = 1;
% time = floor(timeUN)
Inan = find(~isnan(time));  %% added

Igap = find(diff(time(Inan))>time_gap);
Igap(end+1) = length(time(Inan)); %% modified

Igap = Inan(Igap);          %% added

