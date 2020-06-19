function [a,b,c]=planefit(x,y,z)
%function [a,b,c]=bilinreg(x,y,z)
%BILINREG Calculates Bi-Linear-Regrassion
%   [A,B,C]=BILINREG(x,y,z), Where x,y,z are vectors of points
%   values [x(1),y(1),z(1) ;  x(2),y(2),z(2)...x(n),y(n),z(n)].
%   The result are the coefficients to the formula:
%
%              Z = Ax + By + C.
%
%   MCB, NRL, 2014-03-19

mat=[ mean(x.^2)   mean(x.*y)    mean(x)  ;
      mean(x.*y)   mean(y.^2)    mean(y)  ;
      mean(x)       mean(y)        1         ];

vec= [ mean(x.*z)
       mean(y.*z)
       mean(z)    ];

res=(inv(mat)*vec);
a=res(1);
b=res(2);
c=res(3);