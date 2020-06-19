%% function [lonout,latout,isel,jsel] = find_nearest_point_hycom(plon,plat,xpnt,ypnt)
%  Maarten Buijsman, USM, 2016-03-15
%  for a point it finds the nearest grid point on HYCOM grid 
%  input:  plon,plat from HYCOM
%  output: lonout,latout,isel,jsel coordinates and indexes of plon, plat

function [lonout,latout,isel,jsel] = find_nearest_point_hycom(plon,plat,xpnt,ypnt);

[ny,nx] = size(plon);

 % get input by clicking on figure
%[xpnt,ypnt] = ginput(1);

%holder
%plot(xpnt,ypnt,'g*')

% compute distance between point and grid coordinates
% spheric_dist is another function
dist = spheric_dist(plat,ypnt*ones(size(plat)),plon,xpnt*ones(size(plat)));

% find jsel (along y) and isel (along x) index values 
[II,JJ]= meshgrid(1:nx,1:ny);
ii = reshape(II,1,nx*ny);
jj = reshape(JJ,1,nx*ny);
la = reshape(plat,1,nx*ny);
lo = reshape(plon,1,nx*ny);
dd = reshape(dist,1,nx*ny);

% weirdly it finds two minima, on the on northern 
% and one one the southern hemisphere
if ypnt>=0
    Ip = find(la<0);
    ii(Ip)=[];
    jj(Ip)=[];
    la(Ip)=[];
    lo(Ip)=[];    
    dd(Ip)=[];    
elseif ypnt<=0
    Ip = find(la>0);
    ii(Ip)=[];
    jj(Ip)=[];
    la(Ip)=[];
    lo(Ip)=[];    
    dd(Ip)=[];    
end
[mindist,I] = min(dd);

% these values are used in the Fortran diagnostics
isel = ii(I);
jsel = jj(I);

%plot(plon(jsel,isel),plat(jsel,isel),'rs')

% for output to main code
lonout = plon(jsel,isel);
latout = plat(jsel,isel);

% output to screen
disp(' ')
disp(sprintf('jsel = %5i',jsel));
disp(sprintf('isel = %5i',isel));

disp(' ')
disp(sprintf('lon = %8.3f',plon(jsel,isel)));
disp(sprintf('lat = %8.3f',plat(jsel,isel)));



% uncomment this to plot the two minimums
% figure
% pcolor(plon(1:stp:end,1:stp:end),plat(1:stp:end,1:stp:end),dist(1:stp:end,1:stp:end)); shading flat
% colorbar
% hold
% plot(plon(jsel,isel),plat(jsel,isel),'rs')
