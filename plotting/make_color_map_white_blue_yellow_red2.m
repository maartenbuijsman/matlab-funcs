%% make_color_map_white_blue_yellow_red2.m
%% Maarten Buijsman, GFDL, 2011-11-12
%% dark blue blue white yellow red and dark red

function [CMAP] = make_color_map_white_blue_yellow_red2;

jet2 = jet(64);
CMAP1 = jet2; 

[l,w] = size(jet2);
 
%Inum = 25;
Inum = 40;

%% remove bottom part
CMAP1(1:Inum,:)=[];

%% intrep bottom part
CMAP2=[];
white = [1 1 1];
CMAP2(1:floor(Inum/2),:) = interp1([1 floor(Inum/2)],[white; CMAP1(1,:)],1:floor(Inum/2));
%CMAP2(1:Inum,:) = interp1([1 Inum+1],[white; CMAP1(1,:)],1:Inum);

CMAP = [CMAP2; CMAP1];

% %% test
% figure
% colormap(CMAP)
% colorbar