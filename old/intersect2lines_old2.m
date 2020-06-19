%% function [xi,yi,isInSegment] = intersect2lines(x,y);
%% MCB, UCLA, 2010-06-15
%% http://stackoverflow.com/questions/2050850/matlab-find-point-of-intersection-between-two-vectors
%% computes intersection of two lines xi, yi, and returns isInSegment=1 if value is between endpoints
%% input are the x and y begin and endpoints of the lines
%% x = [x1(1) x2(1); x1(end) x2(end)];  %# Starting points in first row, ending points in second row
%% y = [y1(1) y2(1); y1(end) y2(end)];

function [xi,yi,isInSegment] = intersect2lines(x,y);

dx = diff(x);  %# Take the differences down each column
dy = diff(y);
den = dx(1)*dy(2)-dy(1)*dx(2);  %# Precompute the denominator
if den~=0
    ua = (dx(2)*(y(1)-y(3))-dy(2)*(x(1)-x(3)))/den;
    ub = (dx(1)*(y(1)-y(3))-dy(1)*(x(1)-x(3)))/den;
    
    xi = x(1)+ua*dx(1);
    yi = y(1)+ua*dy(1);
else
    ua =  
    xi = x(1)+ua*dx(1);
    yi = y(1)+ua*dy(1);
end


isInSegment = all(([ua ub] >= 0) & ([ua ub] <= 1));

return

%% determine a's and b's of y=ax+b



 I = inv(m_1 -1; m_2 -1]*[-b_1;-b_2];

X = A\B
 

X = A\B
