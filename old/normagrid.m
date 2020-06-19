%%---------------------------------------------------------------------------------------------------
%%-- function [norm_out,grid_out] = normagrid(data,norm_grid);
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 24-02-2004
%%-- Date 1st programmed: 18-08-2003
%%-- Description: 
%%-- normalizes data for gridding purp. in norm_grid_calc.m.
%%-- subtracts minimum and normalizes by amplitude
%%-- vector norm_grid is grid to be gridded on 
%%-- vector data is input data
%%-- vector norm_out is normalized data
%%-- vector grid_out is of same length as norm_grid and can be used for
%%-- contour and imagesc
%%---------------------------------------------------------------------------------------------------
function [norm_out,grid_out] = normagrid(data,norm_grid);

norm_out = (data - min(data))/(max(data)-min(data)); % for gridding
grid_out = linspace(min(data),max(data),length(norm_grid)); % for plotting
