function [xx,yy] = stairsc(xf,xc,yc)

%% function [xx,yy] = stairsc(xf,xc,yc)
%  mcb, USM, 2019-3-15
%  creates stairs plot, where the center yc values (at xc)
%  are mapped to the faces at xf

xx=[]; yy=[];
for i=1:length(yc)
    xx = [xx xf(i) xf(i+1)];
    yy = [yy yc(i) yc(i)];
end