%% function dxd = stretch_grid(numcell,maxpct,numcelltran,dxc)
%% Maarten Buijsman, GFDL, 2011-05-27
%% numcell=number of cells in stretched grid > numcelltran
%% maxpct=max percentage
%% numcelltran=number of cells in transition
%% dxc=undisturbed center grid cell
%% dxd=DX
%% NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE
%% NOTE numcelltran*(1+1/fac) needs to be an integer number

function dxd = stretch_grid(numcell,maxpct,numcelltran,dxc);

%test
% numcell = 100;
% maxpct = 0.04;
% numcelltran = 50;   
% dxc=DX;

%% all dzs indicate percent increase in dx
dzsur = 0;        
%maxpct = 0.04;  %% max percentage
%        maxpct = 0.017;  %% max percentage
%        maxpct = 0.017;  %% max percentage
dzi   = maxpct-dzsur;
sp    = 1; %% scaling paramter to determine speed of 

Hzero = 3;
%numcelltran = 60;      %% number of cells over which transition occurs
nz2  = numcelltran;
z1   = linspace(1,numcelltran,nz2);
dz1  = dzi*(tanh( sp*(z1*(2*pi)/numcelltran-pi) )+1)/2 + dzsur;

%% do some scaling
dz1 = dz1-dz1(1);
dz1 = dz1/max(dz1)*dzi;
ndxside_left = numcell - numcelltran; %% these cells have constant percent increase 

%% and back to zero percent => to get a few constant size grid cells
%% added faster interpolation
fac = 3;
if numcelltran*(1+1/fac)+Hzero < numcell
    %% from inside to wall boundary
    dz2 = [dz1 dz1(end)*ones(1,numcell-Hzero-numcelltran*(1+1/fac)) ...
        interp1(1:numcelltran,dz1(end:-1:1),linspace(1,numcelltran,numcelltran/fac)) zeros(1,Hzero)]; 
else
    disp('gridding error')
    return
end

%figure; plot(dz2,'k.-')

%% compute dx based on % increase, from left to right
dxd=[];
dxs = dxc;
for i=1:numcell
    dxd(i) = dxs + dxs*dz2(i);
    dxs = dxd(i);
end

return