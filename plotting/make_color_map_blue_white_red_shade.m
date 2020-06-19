%% make_color_map_blue_white_red_shade.m
%% Maarten Buijsman, GFDL, 2011-11-12
%% dark blue blue white yellow red and dark red

function [CMAP] = make_color_map_blue_white_red_shade;

yel2 = [1 1 1]; 
yel  = [1 1 1];
R_c = [1 0 0];
B_c = [0    0.0625    1.0000];
B_d = [0         0    0.5312];
C_d = [0.5000    0         0];


num = 32;

d1 = linspace(B_d(1),B_c(1),num);
d2 = linspace(B_d(2),B_c(2),num);
d3 = linspace(B_d(3),B_c(3),num);
r1 = linspace(B_c(1),yel(1),num);
r2 = linspace(B_c(2),yel(2),num);
r3 = linspace(B_c(3),yel(3),num);
b1 = linspace(yel2(1),R_c(1),num);
b2 = linspace(yel2(2),R_c(2),num);
b3 = linspace(yel2(3),R_c(3),num);
c1 = linspace(R_c(1),C_d(1),num);
c2 = linspace(R_c(2),C_d(2),num);
c3 = linspace(R_c(3),C_d(3),num);



C1 = [d1 r1 b1 c1];
C2 = [d2 r2 b2 c2];
C3 = [d3 r3 b3 c3];

CMAP = [C1;C2;C3]';
