%%---------------------------------------------------------------------------------------------------
%% Function: [XX,YY] = coord_trans(X,Y,alfa)
%% Author: Maarten Buijsman
%% Place: NIOZ
%% Date: 25-09-2003
%% Date 1st programmed: 24-06-2003
%% Description: 
%% Transforms coordinates
%% Input: X (vector), Y(vector), and alfa (radii) 
%% Output: transformed XX (vector) and YY(vector)
%%---------------------------------------------------------------------------------------------------
function [XX,YY] = coord_trans(X,Y,alfa)

XX = X.*cos(alfa) + Y.*sin(alfa);
YY = -X.*sin(alfa) + Y.*cos(alfa);
