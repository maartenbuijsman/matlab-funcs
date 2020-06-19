%% function [time_new] = insertNaN(time,Igap);
%% Maarten Buijsman, NIOZ, 26-09-04
%% function inserts NaNs in time gaps 
%% uses function find_gap to define time gaps in vector time series 
%% input vector time and time_gap indicating 1 hour, 1 day, etc
%% makes plots look better
function [time_new] = insertNaN02(time,time_gap);

[Igap] = find_gap(time,time_gap)

istart = 1;
time_new = [];
for i=1:length(Igap)
    time_store = [time_new time(istart:Igap(i))];
    time_new = [time_store NaN];

    istart = Igap(i)+1;
end
