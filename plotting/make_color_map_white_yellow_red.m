%% make_color_map_white_yellow_red.m
%% Maarten Buijsman, GFDL, 2011-11-12
%% dark blue blue white yellow red and dark red

function [CMAP] = make_color_map_white_yellow_red;

% B_d  = [1 1 1]; 
% yel2 = [1 1 0];
% 
% R_c = [1 0 0];
% B_c = [0    0.0625    1.0000];
% C_d = [0.5000    0         0];
% 
% 
% num = 32;
% 
% d1 = linspace(B_d(1),yel2(1),num); % white-yellow
% d2 = linspace(B_d(2),yel2(2),num);
% d3 = linspace(B_d(3),yel2(3),num);
% b1 = linspace(yel2(1),R_c(1),num); % yellow-red
% b2 = linspace(yel2(2),R_c(2),num);
% b3 = linspace(yel2(3),R_c(3),num);
% c1 = linspace(R_c(1),C_d(1),num);  % red-shaded
% c2 = linspace(R_c(2),C_d(2),num);
% c3 = linspace(R_c(3),C_d(3),num);
% 
% 
% C1 = [d1 b1];
% C2 = [d2 b2];
% C3 = [d3 b3];
% 
% CMAP = [C1;C2;C3]';

jet2 = jet;

B_d  = [1 1 1]; 
num  = 24;

jet2(1:num,1) = linspace(B_d(1),jet2(num,1),num)';
jet2(1:num,2) = linspace(B_d(2),jet2(num,2),num)';
jet2(1:num,3) = linspace(B_d(3),jet2(num,3),num)';

CMAP = jet2;