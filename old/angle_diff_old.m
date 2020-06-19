%% function [ang_dif] = angle_diff(ang1,ang2)
%% MCB, UCLA, 19-03-2008
%% compute diffrence between two angles = ang1 - ang2
%% input: vectors of angles in columns ang1 ang2 in radians
%% outtput: difference angle alp in radians
%% HOWEVER, the answer is always positive

function [alp] = angle_diff(ang1,ang2)

zabs = sqrt(2); %% length of vector

%% get x and y of vectors
x1 = zabs*cos(ang1); y1 = zabs*sin(ang1);    
x2 = zabs*cos(ang2); y2 = zabs*sin(ang2);    

alp = acos(dot([x1'; y1'],[x2'; y2'])/(zabs*zabs))'; %% improduct

%% test
% psi*180/pi
% ang1=psi(:,1); ang2=psi(:,2);
% plot([0 x1],[0 y1],'b-',[0 x2],[0 y2],'r-')
% axis equal
% alp*180/pi


