%% function [xx2,yy2,zz2] = grid_half(xx,yy,zz);
%% IGPP, Maarten Buijsman 2007-10-18
%% calculates the half-way positions
%% for plotting with pcolor

function [xx2,yy2,zz2] = grid_half(xx,yy,zz);

%% test
% xx = [1 2 3 4;1 2 3 4;1 2 3 4];
% yy = [-10 -8 -8 -10;-5 -4 -4 -5;0 0 0 0];
% zz = [3.1 3.2 3.3 3.4; 2.1 2.2 2.3 2.4; 1.1 1.2 1.3 1.4];

%% another test --------------------------
% x=[0:10]; y=[0:10];
% [xx,yy] = meshgrid(x,y);
% 
% %zz=ones(size(x))'*[cos(2*pi/10*x)];
% zz=cos(2*pi/10*y)'*[cos(2*pi/10*x)];
% 
% figure
% contourf(xx,yy,zz,40)
% colorbar
% hold
% 
% %plot(x,cos(2*pi/10*x))
% amp=ones(size(x))'*[1-0.75*sin(2*pi/20*x)];
% yyy = yy.*amp;
% zzz=cos(2*pi/10*yyy).*[cos(2*pi/10*xx)];
% yy=yyy;
% zz=zzz;
%% another test --------------------------


%% find the centers of the cells
xx2 = (xx(1:end-1,1:end-1)+xx(2:end,1:end-1)+xx(1:end-1,2:end)+xx(2:end,2:end))/4;
yy2 = (yy(1:end-1,1:end-1)+yy(2:end,1:end-1)+yy(1:end-1,2:end)+yy(2:end,2:end))/4;

%% throw away the first row and column, which can not be mapped anyway.
zz2 = zz(2:end,2:end);

%% another test --------------------------
%figure; pcolor(xx2,yy2,zz2)
% pcolor(xx2,yy2,zz2)
% shading flat
%% another test --------------------------
