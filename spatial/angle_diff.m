%% function [ang_dif] = angle_diff(ang1,ang2)
%% MCB, UCLA, 15-04-2008
%% compute difference between two angles = ang1 - ang2
%% 
%% input: vectors of angles in columns ang1 ang2 in radians
%% outtput: difference angle alp in radians
%% HOWEVER, the answer is always positive

function [alp] = angle_diff(ang1,ang2);

% clear i
% ang1=-180*pi/180
% ang2=10*pi/180

z = exp(i*(ang1-ang2));
alp = atan2(imag(z), real(z));           %radians
%alp = atan2(imag(z), real(z))*180/pi;   %degrees

