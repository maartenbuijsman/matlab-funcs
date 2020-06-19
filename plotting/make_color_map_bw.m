%%---------------------------------------------------------------------------------------------------
%%-- Program: make_color_map_bw.m
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 30-3-2004
%%-- Description: 
%%-- Makes black and white colormaps

function [CMAP] = make_color_map_bw;

%% make blue red colormap
%yel2 = [0.99 0.99 1/1.5]; 
%yel = [0.99 0.99 0.99];   
%yel2 = [245 245 190]/255;
%yel  = [225 240 240]/255;

yel  = [1 1 1];
B_c = [0 0 0];

num = 64;
r1 = linspace(B_c(1),yel(1),num);
r2 = linspace(B_c(2),yel(2),num);
r3 = linspace(B_c(3),yel(3),num);

C1 = [r1];
C2 = [r2];
C3 = [r3];

%% trick
C1 = [1 0.8]; C2 = [1 0.8]; C3 = [1 0.8];
CMAP = [C1;C2;C3]';
