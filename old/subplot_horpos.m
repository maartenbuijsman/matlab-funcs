%% function [pos] = subplot_horpos(numplot,hors,vers,Ds)
%% Maarten Buijsman, NIOZ, 08-10-04
%% determines position horizontal subplots

function [pos] = subplot_horpos(numplot,hors,vers,Ds);
Dv = 1-2*vers;
Dh = (1-(numplot)*Ds-hors)/numplot;
pos(1,:)=[hors vers Dh Dv];
for i=2:numplot
    pos(i,:)=[hors+(Dh+Ds)*(i-1) vers Dh Dv];
end
