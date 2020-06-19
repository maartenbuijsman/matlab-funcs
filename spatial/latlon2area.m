%% function A = latlon2area(lon,lat);
%  USM, MCB, 2020-06-2020
%  compute area given vector lon and lat
%  of cartesian coordinates
%  area is computed at coordinates
%  in funcs/old/area_grid.m the geoid is used, but
%  area calculation has a singularity

function A = latlon2area(lon,lat);

%lon = lo;lat = la; %test

[LO,LA] = meshgrid(lon,lat);
[ny,nx] = size(LO);

if abs(mean(diff(LA(:,1)))-diff(LA(1:2,1)))>1e-10;
    mean(diff(LA(:,1)))
    diff(LA(1:2,1))
    disp('ERROR; not equidistant');
    pause
end

% get area based on xf and yf
dx = spheric_dist(LA(:,1),LA(:,2),LO(:,1),LO(:,2));
dy = spheric_dist(LA(1,1),LA(2,1),LO(1,1),LO(2,1));
A  = repmat(dx*dy,[1 nx]);
