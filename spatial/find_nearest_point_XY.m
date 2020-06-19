function [Xg,Yg,Ix,Iy] = ffind_nearest_point_XY(Xw,Yw,Xgin,Ygin);

%% [Xg,Yg,Ix,Iy] = find_nearest_point_XY(Xw,Yw,Xgin,Ygin);
% USM, MCB, 2020/05/02
% find nearest point on a cartesian grid (in x,y coordinates)
% input are grid matrices Xw and Yw 
% Xgin,Ygin are points that should be matched to grid locations
% formerly known as find_nearest_points_2020_v1.m

%test
%Xw=X; Yw=Y; Xgin=xp; Ygin=yp;
  
%% find nearest grid point

for i=1:length(Xgin)

    % compute distances 
    dx = Xgin(i) - Xw;
    dy = Ygin(i) - Yw;
    ds = sqrt(dx.^2 + dy.^2);
    
    % find indices
    mindistx = min(ds,[],1); % along x
    [mindistx2,Ix(i)] = min(mindistx);
    
    mindisty = min(ds,[],2); % along y
    [mindisty2,Iy(i)] = min(mindisty);    
    
    % map indices to coordinates
    Xg(i) = Xw(Iy(i),Ix(i));
    Yg(i) = Yw(Iy(i),Ix(i));
    
end

