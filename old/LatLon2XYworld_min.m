%% function [Xw,Yw] = LatLon2XYworld_min(lat,lon);
%% Maarten Buijsman, NIOZ, 21-07-06, UCLA, 17-03-2008
%% calculate lat lon in meters from 0 north and 0 west
%% from http://www.ncgia.ucsb.edu/education/curricula/giscc/units/u014/u014.html)
%% computes the distances relative to left and lower corner of (lon,lat)

function [Xw,Yw] = LatLon2XYworld_min(lat,lon);

% Length of a degree of longitude = cos (latitude) * 111.325 kilometers

dist = 111325; %meters
Xw = (lon-min(lon)).*cos(lat*pi/180)*dist;
Yw = (lat-min(lat))*dist;
