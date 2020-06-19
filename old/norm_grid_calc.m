%%---------------------------------------------------------------------------------------------------
%%-- function [Xgrid,Ygrid,Zgrid] = norm_grid_calc(Xin,Yin,Zin,num_grid);
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 25-02-2004
%%-- Date 1st programmed: 25-02-2004
%%-- Description: 
%%-- normalizes data and grids data
%%-- input are vectors Xin, Yin, and Zinof time-series with equal length
%%-- num_grid = number of horizontal and vertical grid points
%%-- Xgrid & Ygrid are vectors for plotting
%%-- Zgrid is the gridded data matrix
%%-- use in contour(Xgrid,Ygrid,Zgrid,[..],'k')
%%-- str prescribes 'linear' or 'cubic'
%%---------------------------------------------------------------------------------------------------
function [Xgrid,Ygrid,Zgrid] = norm_grid_calc(Xin,Yin,Zin,num_grid,str);

%% grid data
%% normalize grid data between 0 and 1
norm_grid = linspace(0,1,num_grid);
[NX,NY] = meshgrid(norm_grid,norm_grid);
[Xnorm,Xgrid] = normagrid(Xin,norm_grid);
[Ynorm,Ygrid] = normagrid(Yin,norm_grid);
Zgrid = griddata(Xnorm,Ynorm,Zin,NX,NY,str);
