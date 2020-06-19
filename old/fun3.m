function out=fun3(coeff,X,Y)
a = coeff(1);
b = coeff(2);
c = coeff(3);
d = coeff(4);
e = coeff(5);
aa = coeff(6);
cc = coeff(7);
dd = coeff(8);

bb = 0;
ff = 0;
gg = 0;
hh = 0;

%Y_fun = (dd + a*X + e*X.^2 + aa*X.^3)...
%      ./(cc + b*X + c*X.^2 +  d*X.^3);
  

%bb = coeff(9); ff = coeff(10); gg = coeff(11); hh = coeff(12);

Y_fun = (dd + a*X + e*X.^2 + aa*X.^3 + bb*X.^4 + gg*X.^5)...
       ./(cc + b*X + c*X.^2 +  d*X.^3 + ff*X.^4 + hh*X.^5);

DIFF = Y_fun - Y;
SQ_DIFF = DIFF.^2;

out = sum(SQ_DIFF);
