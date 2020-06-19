%% function [time_new] = insertNaN(time,Igap);
%% Maarten Buijsman, NIOZ, 5-2-04
%% function inserts NaNs in time gaps defined with
%% function find_gap in vector series time
%% makes plots look better
%% used in combination with funct. find_gap
function [time_new] = insertNaN(time,Igap);

%% make sure it is a vector
[vl,hl] = size(time);
if vl>hl; time = time'; end

istart = 1;
time_new = [];
for i=1:length(Igap)
    time_store = [time_new time(istart:Igap(i))];
    time_new = [time_store NaN];

    istart = Igap(i)+1;
end

%% stuff after Igap
time_new = [time_new time(Igap(i)+1:end)];
