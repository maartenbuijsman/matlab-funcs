%% function [pos] = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);
%% Maarten Buijsman, NIOZ, 06-10-05
%% determines positions of subplots
%% numph = number of horizontal plots, numph = number of vertical plots,
%% hors = start spacing horizontal, hore = end spacing horizontal
%% vers = start spacing vertical, vere = end spacing vertical
%% Dsh = horizontal distance between the plots, Dsv = vertical distance between the plots
%% pos follows the numbering of subplots (left to right, top to bottom) 
%% 
%% example
%% numph = 3;numpv = 3;
%% hors = 0.1;hore =0.1
%% vers = 0.1;vere =0.1
%% Dsh = 0.05;Dsv =0.05
%% [pos] = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);
%%
%% use in combination with 
%% f1 = figure
%% set(f1,'PaperUnits','centimeters','PaperPosition',[0.5 0.5 20 15]);


function [pos] = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);

%% Dh,Dv per plot
Dh = (1-(numph-1)*Dsh-hors-hore)/numph;
Dv = (1-(numpv-1)*Dsv-vers-vere)/numpv;
pos = []; hh=0;
for j=1:numpv
    for i=1:numph
        hh = hh+1;
        pos(hh,:)=[hors+(i-1)*(Dh+Dsh) 1-(vere+Dv+(Dv+Dsv)*(j-1)) Dh Dv];
%         subplot('position',pos(hh,:))
    end
end

