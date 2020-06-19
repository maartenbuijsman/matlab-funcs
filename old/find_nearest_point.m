%% [Xg,Yg,Ix,Iy] = find_nearest_point(Xw,Yw,Xgin,Ygin,str);
%% UCLA, MCB, 2008/03/14
%% find nearest point (in x,y coordinates)
%% input are grid matrices Xw and Yw 
%% Xgin,Ygin are points that should be matched to grid locations
%% input str='map' for finding points in a topographic map
%% else use something else 'nomap'

function [Xg,Yg,Ix,Iy] = find_nearest_point(Xw,Yw,Xgin,Ygin,str);

%% correct for scaling factor
if strcmp(str,'map')
    if     round(floor(mean(mean(Xw)))/Xgin/1e3)==1; fac = 1e3;
    elseif round(floor(mean(mean(Xw)))/Xgin*1e3)==1; fac = 1e-3;
    else fac=1;
    end
else fac=1;
end

Xw = Xw/fac;
Yw = Yw/fac;
   
%% find nearest grid point
ds_min=ones(1,size(Yw,2)); 
IIy=ds_min; 
for ii=1:size(Yw,2);
    dx = abs(Xw(:,ii)-Xgin);
    dy = abs(Yw(:,ii)-Ygin);    
    ds = sqrt(dx.^2+dy.^2);
    [ds_min(ii),IIy(ii)] = min(ds); %min distance
%    plot(Xw(IIy(ii),ii),Yw(IIy(ii),ii),'g.')
end

%whos IIy ds_min
%figure; plot(Xw(IIy,:),ds_min,'k+',Xgin,0,'ro')

%% locate minimum distance
[d,Imin] = min(ds_min); %Imin is along x
%plot(Xw(IIy(Imin),Imin),Yw(IIy(Imin),Imin),'ks')

%% for output
Ix = Imin;
Iy = IIy(Imin);
Xg = Xw(IIy(Imin),Imin);
Yg = Yw(IIy(Imin),Imin);