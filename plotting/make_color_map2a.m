%%---------------------------------------------------------------------------------------------------
%%-- Program: make_color_map2a.m
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 30-3-2004
%%-- Description: 
%%-- Makes colormaps
%% thi sone is symmetrical (2 is not)

function [CMAP] = make_color_map2a;

%% make blue red colormap
%yel2 = [0.99 0.99 1/1.5]; 
%yel = [0.99 0.99 0.99];   
%yel2 = [245 245 190]/255;
%yel  = [225 240 240]/255;

yel2 = [255 255 255/1.5]/255;
yel  = [1 1 1];
R_c = [1 0 0];
B_c = [0    0.0625    1.0000];

num = 64;
r1 = linspace(B_c(1),yel(1),num);
r2 = linspace(B_c(2),yel(2),num);
r3 = linspace(B_c(3),yel(3),num);
b1 = linspace(yel2(1),R_c(1),num);
b2 = linspace(yel2(2),R_c(2),num);
b3 = linspace(yel2(3),R_c(3),num);

C1 = [r1 b1];
C2 = [r2 b2];
C3 = [r3 b3];

CMAP = [C1;C2;C3]';
