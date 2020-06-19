 function out=fun2(coeff,X,Y)
 a = coeff(1);
 b = coeff(2);
 c = coeff(3);
 Y_fun = a .* X + b .* sin(X)+c;
 DIFF = Y_fun - Y; 
 SQ_DIFF = DIFF.^2;

 out = sum(SQ_DIFF);
