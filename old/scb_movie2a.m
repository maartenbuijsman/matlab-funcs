%---------------------------------------------
% make movies of cross-sectional T/S field 
% for SCB 1 km result
%---------------------------------------------

close all;
clear all;

dampmovie=true;
plot_sec=true;

casen=2; % 1:off SMB, 2:
romsgrid='scb_iswake_grid_2.nc';
coastline='/home/yusuke/matlab/Roms_tools/My_tools/usw_coast_h.mat';
mov_name=['/scratch/yusuke/roms/tmp/scb_cross_t' int2str(casen) '.avi'];
compress=['/scratch/yusuke/roms/tmp/scb_cross_t' int2str(casen) '.mp4'];
transfig=['figures/transect' int2str(casen) '.png'];
hisdir='case0/';
fps=8; fr=0;
vrate='1000000';	% video bitrate (unit: bps)
vscale='640:640';	% video aspect ratio (multiple of 8)
vcodec='msmpeg4v2';	% video codec: use msmpeg4v2 for windows xp compatibility

cax=[11.5 19.5]; kvar=40; Defs=12;
bright =load('/home/yusuke/etc/ncview/bright.ncmap')/255;
rainbow=load('/home/yusuke/etc/ncview/rainbow.ncmap')/255;
cmap=rainbow;

dc=(cax(2)-cax(1))/128.0;

nc=netcdf(romsgrid,'nowrite');
lonr=nc{'lon_rho'}(:); 
latr=nc{'lat_rho'}(:); 
h=nc{'h'}(:); 
mask=nc{'mask_rho'}(:); 
close(nc);
mask(find(mask==0))=NaN;

corr=find(lonr<0); 
lonr(corr)=lonr(corr)+360.0;
lonmin=min(min(lonr)); lonmax=max(max(lonr));
latmin=min(min(latr)); latmax=max(max(latr));

%lonsec1=-118.6; lonsec2=-119.8; latsec1=34.04; latsec2=32.7;
if (casen==1);
lonsec1=-119.8; lonsec2=-118.6; latsec1=32.7; latsec2=34.04; dist=185; dmin=-1400;
elseif (casen==2);
lonsec1=-118.4; lonsec2=-117.7; latsec1=32.55; latsec2=33.47; dist=120; dmin=-1200;
end;
lonsec=[lonsec1 lonsec2]; latsec=[latsec1 latsec2];

if (plot_sec);
cax=[0 4000];
hF=figure('Units','centimeters','Position',[5 5 17.5 13.5],...
   'Color','w','PaperPositionMode','auto');
axes('Units','centimeters','Position',[1.5 1 15 0.5]);
cbx=[0:1]; cby=[cax(1):dc:cax(2)]; [CBX,CBY]=meshgrid(cbx,cby);
pcolor(CBY,CBX,CBY); shading('flat'); colormap(cmap);
set(gca,'YTick',cby,'YTickLabel',[' '],'fontsize',Defs);
axes('Units','centimeters','Position',[1.5 1.5 15 12]);
m_proj('Equidistant Cylindrical','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(lonr,latr,h); shading('flat'); caxis(cax); hold on;
m_usercoast(coastline,'patch',[.7 .7 .7],'edgecolor','k'); hold on;
m_grid('box','fancy','tickdir','out','fontsize',Defs);
hL=m_line(lonsec+360,latsec); set(hL,'linewidth',2,'color','w');
title('bathymetry (m)');
print('-dpng','-r150',transfig);
end;

tic;
cax=[5 20];
if (dampmovie); ifr=1:90; trg=1:8; else; ifr=90:90; trg=1:1; end;

for ifile=ifr;

ncfile=[hisdir 'scb_his.' my_get_date((ifile-1)*8,4) '.nc'];
disp(['processing ' ncfile]);

for tind=trg;

fr=fr+1; count=(ifile-1)*8+tind-1;
iday=ceil(count/8)+1; ihour=mod(count,8)*3;
stime=[my_get_date(iday,3) ' day ' my_get_date(ihour,2) ' hour (2002)'];
[x,z,dat]=get_section(ncfile,romsgrid,lonsec,latsec,'temp',tind);

hF=figure('Units','centimeters','Position',[5 5 17.5 17.5],...
   'Color','w','PaperPositionMode','auto');

axes('Units','centimeters','Position',[2 8.5 15 8]);
pcolor(x,z,dat); shading('flat'); caxis(cax); hold on;
axis([0 dist -200 2]);
title(['temperature (deg. C) : ' stime]); ylabel('depth (m)'); 
set(gca,'xticklabel',[],'ytick',[-180:20:0]);

axes('Units','centimeters','Position',[2 1.5 15 6.9]);
pcolor(x,z,dat); shading('flat'); caxis(cax); hold on;
axis([0 dist dmin -200]);
fill([0 0 0 0],[0 0 0 0],'k');
xlabel('along-transect distance (km)'); ylabel('depth (m)'); 

axes('Units','centimeters','Position',[9.5 2 7 0.5]);
cbx=[0:1]; cby=[cax(1):dc:cax(2)];
[CBX,CBY]=meshgrid(cbx,cby);
pcolor(CBY,CBX,CBY); shading('interp'); colormap(cmap);
set(gca,'YTick',cby,'YTickLabel',[' '],'fontsize',Defs,...
    'tickdir','out','xaxislocation','top');

if (dampmovie); M(fr)=getframe(hF); close all; end;

end;    % end loop on tind

end;    % end loop on ifile

if (dampmovie); 
movie2avi(M,mov_name,'Fps',fps); 
system(['mencoder ' mov_name ' -ovc lavc -lavcopts vcodec=' vcodec ':vbitrate=' vrate ...
        ' -vf scale=' vscale ' -o ' compress]);
system(['/bin/rm ' mov_name]);
end;
toc;

return;

