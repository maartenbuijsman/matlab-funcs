%% function [x_s,y_fit,r,p,cf] = line_fit(x,y,N);
%% Maarten Buijsman, NIOZ, 10-06-05
%% fits lines through data, and produces correlation coefficient
%% input: vectors x, y, and the order N
%% output: vectors x_s, y_fit, corrcoef r, and coefficients P
function [x_s,y_fit,r,p,cf] = line_fit(x,y,N);

Inan = find(isnan(x)); x(Inan) = []; y(Inan) = [];
Inan = find(isnan(y)); x(Inan) = []; y(Inan) = [];

%% sort
[x_s,Is]=sort(x);
y_s = y(Is);

% %% normalize x_s
% x_n = (x_s - mean(x_s))./std(x_s);
% [cf] = polyfit(x_n,y_s,N);
% [y_fit] = polyval(cf,x_n);
% % %% correct coefficients %% does not work ....
% cf(1) = cf(1)*std(x_s);
% cf(2) = cf(2)-cf(2)*mean(x_s);

%% not normalized produces correct coefficients a and b of the line
[cf] = polyfit(x_s,y_s,N);
[y_fit] = polyval(cf,x_s);

%% correlation coefficient
[RR,PP] = corrcoef(x,y); 
p = PP(2,1);
r = RR(2,1); 


