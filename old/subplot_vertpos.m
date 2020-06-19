%% function [pos] = subplot_vertpos(numplot,hors,hore,vers,Ds);
%% Maarten Buijsman, NIOZ, 08-10-04
%% determines position vertical subplots
%% numplot = number of plots, hors = start spacing horizontal
%% hore = end spacing horizontal, vers = start spacing vertical
%% vere = end spacing vertical
%% Ds = distance between the plots

function [pos] = subplot_vertpos(numplot,hors,hore,vers,vere,Ds);
Dh = 1-hors-hore;
Dv = (1-(numplot)*Ds-vers-vere)/numplot;
for i=numplot:-1:2
    pos(numplot-i+1,:)=[hors vers+(Dv+Ds)*(i-1) Dh Dv];
end
pos(numplot,:)=[hors vers Dh Dv];
