%% make_color_map_jet_white.m
%% Maarten Buijsman, GFDL, 2011-11-12
%% dark blue blue white yellow red and dark red

function [CMAP] = make_color_map_jet_white;

jet2 = jet;
CMAP = jet2; 

[l,w] = size(jet2);

Inum = 20;
i1   = Inum;
i2   = l-Inum+1;

jet3 = jet2(1:Inum,:); %% bottom
jet4 = jet2(i2:end,:); %% top

%% now interpolate white
white = [1 1 1];
CMAP(i1:l/2,:) = interp1([i1 l/2],[jet3(end,:); white],i1:l/2);
CMAP(l/2+1:i2,:) = interp1([l/2+1 i2],[white; jet4(1,:)],l/2+1:i2);

% %% test
% figure
% colormap(CMAP)
% colorbar