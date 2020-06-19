%% function [xn,zn] = beam_coords(xs,zs,slope,xdir,zdir,RC,XC,Depth);
%% MCB, GFDL, 2011-08-24
%% computes beam coordinates based on starting point
%%     ixs = 100; start coordinates
%%     izs = 100;
%%     xdir = 'right';  directions
%%     zdir = 'up';
%%     slope = slope1; beam slope

function [xn,zn] = beam_coords(xs,zs,slope,xdir,zdir,RC,XC,Depth);

nx = length(XC);

%% find starting indexes
[d,ixs]=min(abs(XC-xs));
xn = XC(ixs);

[d,izs]=min(abs(RC-zs));
zn = RC(izs);

k=0;
while ixs<nx & ixs>1
    k=k+1;
    if k==1
        if strcmp(xdir,'right'); xcnt=+1;
        else                     xcnt=-1;
        end
        
        if strcmp(zdir,'up'); zcnt=+1;
        else                  zcnt=-1;
        end
    end
    
    ixs = ixs + 1*xcnt;
    xn(k+1) = XC(ixs);
    DX = xn(k+1)-xn(k);
    slopei = interp1(RC,slope,zn(k),'linear','extrap');
    
    DZ = xcnt*zcnt*slopei*DX;
    zn(k+1) = zn(k)+DZ;
    
    if     zn(k+1) > 0;           zn(k+1) = zn(k); zcnt=-zcnt; %flip direction
    elseif zn(k+1) < -Depth(ixs); zn(k+1) = zn(k); zcnt=-zcnt;
    end
    
    
end
