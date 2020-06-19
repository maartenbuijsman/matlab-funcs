%% function [xl,yl,sl,varstruct] = get_slice_mitgcm(XC,YC,varstruct,xlp,ylp,nump,getpoints,scfac);
%% MCB, GFDL, 2012-03-19
%% make a slice of 2D or 3D coordinates
%% given line ends xlp and ylp, and 2D or 3D variables in varstruct
%% slices are make per layer
%% when the points are extracted from figure first do xlp=0; ylp=0; nump=100; getpoints=1;
%% multiply line coordinates by scfac when figure coordinates are in km instead of meters 

function [xl,yl,sl,varstruct] = get_slice_mitgcm(XC,YC,varstruct,xlp,ylp,nump,getpoints,scfac);

%% get points in figure ===============================================
if getpoints
    xlp=[];ylp=[];
    for i=1:2
        [xlp(i),ylp(i)]=ginput(1);
        holder
        plot(xlp(i),ylp(i),'m+')
        drawnow
    end
    plot(xlp,ylp,'k--')
    
    disp(['xlp = [' sprintf('%10.6f ',xlp) '];']);
    disp(['ylp = [' sprintf('%10.6f ',ylp) '];']);
    
end

%% rescale if needed
xlp = xlp*scfac;
ylp = ylp*scfac;

%% get line ===============================================

%% line test
% km
% xlp = [106.031010 118.247671 ]*1e3;
% ylp = [ 38.202536  43.924264 ]*1e3;
% lat,lon
% xlp = [121.927400 122.028043 ];
% ylp = [ 20.638412  20.675561 ];

[minx,Imin] = min(xlp); miny=(ylp(Imin));
[maxx,Imax] = max(xlp); maxy=(ylp(Imax));

%nump=100;
xl=linspace(minx,maxx,nump);
yl=linspace(miny,maxy,nump);
sl=sqrt((xl - xl(1)).^2+(yl - yl(1)).^2);

%% limit data ===============================================
DOFF = XC(1,2)-XC(1,1);  %% max grid size
Ix = find(XC(1,:)>minx-DOFF & XC(1,:)<maxx+DOFF );
Iy = find(YC(:,1)>miny-DOFF & YC(:,1)<maxy+DOFF );

%% grid ===============================================
[tri,ts,ds] = func_init_griddata(XC(Iy,Ix),YC(Iy,Ix),xl,yl);

for ii=1:length(varstruct)
    nz2 = size(varstruct(ii).var,3);
    varstruct(ii).varslice = zeros(nump,nz2);
    
    %% loop per variable
    for i=1:nz2
        disp(['layer ' num2str(i) ' of variable ' num2str(ii)])
        
        varNaN = varstruct(ii).var(Iy,Ix,i);
        if ~isempty(rem_NaN(varNaN(:)))
            varstruct(ii).varslice(:,i) = func_fast_griddata(XC(Iy,Ix),YC(Iy,Ix),varstruct(ii).var(Iy,Ix,i),xl,yl,'linear',tri,ts);
        else
            disp(['    layer ' num2str(i) ' has NaNs only'])            
        end
    end
end

disp(['xlp = [' sprintf('%10.6f ',xlp) '];']);
disp(['ylp = [' sprintf('%10.6f ',ylp) '];']);


return