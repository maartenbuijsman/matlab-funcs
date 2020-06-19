clear all

 x=1:10;
 y=sin(x);
 bestcoeffs=fminsearch(@fun2,[1 1 1],[],x,y);
 yfit=bestcoeffs(1)*x +bestcoeffs(2)*sin(x) +  bestcoeffs(3);
 %Now compare y with yfit

 figure
 plot(x,y,x,yfit);

 
%  function out=fun(coeff,X,Y)
%  a = coeff(1);
%  b = coeff(2);
%  c = coeff(3);
%  Y_fun = a .* X + b .* sin(X)+c;
%  DIFF = Y_fun - Y; 
%  SQ_DIFF = DIFF.^2;
% 
%  out = sum(SQ_DIFF);
%  
%  return

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% x=1:10;
%  y=sin(x);
%  [P,S] = polyfit(x,y,5); 
%  yfit= polyval(P,x);
%  %Now compare y with yfit
% figure
%  plot(x,y,x,yfit);
% 
% %%---------------------------
% 
% function fitting
%  x=1:10;
%  y=sin(x);
%  
%  bestcoeffs=fminsearch(@(x) myfun(x,1,1,1),x);
%  yfit=bestcoeffs(1)*x +bestcoeffs(2)*sin(x) +  bestcoeffs(3);
%  %Now compare y with yfit
%  plot(x,y,x,yfit);
% 
% 
%  
% %%---------------------------
% function fitting
%  x=1:10;
%  y=sin(x);
%  
% myfun= a*x+b*sin(x)+c
%  
%  bestcoeffs=fminsearch(@fun,[1 1 1],[],x,y);
%  yfit=bestcoeffs(1)*x +bestcoeffs(2)*sin(x) +  bestcoeffs(3);
%  %Now compare y with yfit
%  plot(x,y,x,yfit);
% 
%  
%  
% banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
% [x,fval] = fminsearch(banana,[-1.2, 1])
% 
% 
% 
% fminsearch(@cos,[2])
% 
% myfun([0.9298    0.0312]*1.0e-03,1.5)
% 
% myfun = @(x,y,a) x^2 + a*y^2;
% 
% a = 1.5; % define parameter first
% x=0; y=1;
% x = fminsearch(@(x,y) myfun(x,y,a),0)
% 
% function fitting = myfun(x,a)
% fitting = x(1)^2 + a*x(2)^2;
% 
% sqr = @(x) x.^2;
% 
% sqr(15)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  x=1:10;
%  y=sin(x);
%  
%  myfun= @(x,a) a(1)*x+a(2)*sin(x)+a(3)
%  
%  bestcoeffs=fminsearch(@(x,a) myfun(x,a),[1 1 1],[],x,y);
%  yfit=bestcoeffs(1)*x +bestcoeffs(2)*sin(x) +  bestcoeffs(3);
%  %Now compare y with yfit
%  plot(x,y,x,yfit);
% 
