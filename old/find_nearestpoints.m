%% [Xg,Yg,Ix,Iy] = find_nearestpoints(Xw,Yw,Xgin,Ygin,nump);
%% UCLA, MCB, 2011/10/03
%% find nearest nump points (in x,y coordinates)
%% input are grid matrices Xw and Yw 
%% Xgin,Ygin is coordinate for which the nearest point should be found
%% IT ONLY WORKS FOR cartesian coords X,Y
$$
function [Xg,Yg,Ix,Iy] = find_nearestpoints(Xw,Yw,Xgin,Ygin,nump);

% %% test %% test %% test %% test 
% Xw = XG2;
% Yw = YG2;
% Xgin=XX(1,1)
% Ygin=YY(1,1)
% 
% figure
% plot(Xw,Yw,'k.')
% hold
% plot(Xgin,Ygin,'rs')
% 
% %% test %% test %% test %% test 
   
%% find nearest grid point
dx = abs(Xw-Xgin);
dy = abs(Yw-Ygin);    
ds = sqrt(dx.^2+dy.^2);

IXM = repmat(1:size(Xw,2),size(Xw,1));
IYM = repmat(1:size(Xw,1),size(Xw,1));

for II=1:nump
    for ii=1:size(Yw,2); %% along x
        [ds_min(ii),IIy(ii)] = min(ds(:,ii)); %min distance along y
    end
    [d,P(II).Ixmin] = min(ds_min); %Imin is along x
    P(II).Iymin     = IIy(P(II).Ixmin);
    
    %% omit in next round
    ds(P(II).Iymin,P(II).Ixmin) = max(ds(:));
    plot()
end
    
    

ds_min=[]; IIy=[]; 
for ii=1:size(Yw,2); %% along x
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