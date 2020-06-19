%% function [dhdx,dhdy,dhdxi,dhdeta] = xy_grad_rho(h,pm,pn,ang2);
%% Maarten Buijsman, UCLA, 2009-04-15
%% computes gradients of h at rho points of ROMS grid 
%% if h is scalar, grad h is vector
%% if h is vector, grad dot h is scalar
%% along xi and eta axes, and along x and y axes
%% dhdx and dhdy are rotated using angle ang

function [dhdx,dhdy,dhdxi,dhdeta] = xy_grad_rho(h,pm,pn,ang2);

DXI  = 1./pm;
DETA = 1./pn;

%along x
DXI1   = (DXI(:,2:end)+DXI(:,1:end-1))/2; %mean DX
dhdxi1 = diff(h,1,2)./DXI1;            %at u points
dhdxi  = ones(size(pm))*NaN;
dhdxi(:,2:end-1) = (dhdxi1(:,2:end)+dhdxi1(:,1:end-1))/2; %mean at rho points
dhdxi(:,1)       = dhdxi1(:,1);            % fill in boundaries
dhdxi(:,end)     = dhdxi1(:,end);

%along y
DETA1   = (DETA(2:end,:)+DETA(1:end-1,:))/2;
dhdeta1 = diff(h,1,1)./DETA1;
dhdeta  = ones(size(pm))*NaN;
dhdeta(2:end-1,:) = (dhdeta1(2:end,:)+dhdeta1(1:end-1,:))/2;
dhdeta(1,:)       = dhdeta1(1,:);            % fill in boundaries
dhdeta(end,:)     = dhdeta1(end,:);

%% rotate relative to x,y coords
[dhdx,dhdy] = coord_trans(dhdxi,dhdeta,-ang2); 
