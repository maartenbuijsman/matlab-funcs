%% function [xi,yi,isInSegment] = intersect2lines(xx1,yy1,xx2,yy2);
%% MCB, UCLA, 2010-06-15
%% http://stackoverflow.com/questions/2050850/matlab-find-point-of-intersection-between-two-vectors
%% computes intersection of two lines xi, yi, and returns isInSegment=1 if value is between endpoints
%% input are four vectors, each vector has 2 x or y values of each line
%% program determines determine a's and b's of line y=ax+b
%% then solves for intersection and checks if intersection is on the line segments
%% check comments page 247
function [xi,yi,isInSegment] = intersect2lines(xx1,yy1,xx2,yy2);

%% test
%xx1=xps; yy1=yps; xx2=[x1([1 end])]; yy2=[y1([1 end])];

%% determine a's and b's of y=ax+b
if diff(xx1)==0; xx1(1)=xx1(1)+1; end %% trick to bypass infinite slope problem
[d,Imx1]=max(xx1); [d,Imn1]=min(xx1);   %% make sure x increases
a1 = (yy1(Imx1)-yy1(Imn1))/(xx1(Imx1)-xx1(Imn1));
b1 = yy1(Imx1) - a1*xx1(Imx1);

if diff(xx2)==0; xx2(1)=xx2(1)+1; end %% trick to bypass infinite slope problem
[d,Imx2]=max(xx2); [d,Imn2]=min(xx2);
a2 = (yy2(Imx2)-yy2(Imn2))/(xx2(Imx2)-xx2(Imn2));
b2 = yy2(Imx2) - a2*xx2(Imx2);

%% matrices a x - y = -b
A = [a1 -1; a2 -1];
B = [-b1 ; -b2];

X  = A\B;
xi = X(1);
yi = X(2);

%% check if points are on the segment
isInSegment = all(xi>=xx1(Imn1) & xi<=xx1(Imx1) & xi>=xx2(Imn2) & xi<=xx2(Imx2));

return
 


