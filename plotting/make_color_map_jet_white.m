%% make_color_map_jet_white.m
%% Maarten Buijsman, GFDL, 2011-11-12
%% dark blue blue white yellow red and dark red

function [CMAP] = make_color_map_jet_white;

jet2 = jet(64);
CMAP = jet2; 

[l,w] = size(jet2);

Inum = 20;
%Inum = 25;
c1   = Inum;
c2   = l-Inum+1;

jet3 = jet2(1:Inum,:); %% bottom
jet4 = jet2(c2:end,:); %% top

%% now interpolate white
white = [1 1 1];
CMAP(c1:l/2,:) = interp1([c1 l/2],[jet3(end,:); white],c1:l/2);
CMAP(l/2+1:c2,:) = interp1([l/2+1 c2],[white; jet4(1,:)],l/2+1:c2);

% %% test
% figure
% colormap(CMAP)
% colorbar