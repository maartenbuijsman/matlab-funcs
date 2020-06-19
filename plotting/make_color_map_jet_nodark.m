%% make_color_map_jet_white_nodark.m
%% Maarten Buijsman, GFDL, 2011-11-12
%% blue white yellow red

function [CMAP] = make_color_map_jet_white_nodark;

jet1 = jet(128);
jet2 = jet1(11:128-10,:);
CMAP = jet2; 


% %% test
% figure
% colormap(CMAP)
% colorbar