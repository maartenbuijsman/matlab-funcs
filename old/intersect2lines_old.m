%% function [p_int] = intersect2lines(xd,yd,yR);
%% MCB, UCLA, 2008-04-15
%% computes intersection of two lines
%% all x and y values are on same grid

function [p_int] = intersect2lines(xd,yd,yR);

m1 = diff(yd)/diff(xd);  % slope of line connecting (xd(1),yd(1)), (xd(2),yd(2))
b1 = yd(1)-xd(1)*m1;    

m2 = diff(yR)/diff(xd);  % slope of line connecting (xR(1),yR(1)), (xR(2),yR(2))
b2 = yR(1)-xd(1)*m2;    


% find intersection of two lines
p_int = [m1 -1; m2 -1]\[-b1; -b2];
