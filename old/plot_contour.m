%%---------------------------------------------------------------------------------------------------
%%-- Fucntion: plot_contour.m
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 10-1-2003
%%-- Description: 
%%-- plots data with contours 
%%-- note that the x series have to be corretced by +x_min 
%%-- and that x series are scaled to produce better xyz
%%-- input: [scaledata,range,x_grid,y_grid,xyz]
%%-- output: - 
%%---------------------------------------------------------------------------------------------------
function plot_contour(scaledata,x_min,range,fontsize,LineWidth,X,Y,xyz);

[C,h] = contour(X/scaledata+x_min,Y,xyz,range,'k');
set(h,'LineWidth',LineWidth);

if length(C)>0
    h2 = clabel(C,h);
    set(h2,'FontSize',fontsize);
end

