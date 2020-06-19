%% function dxd = stretch_grid_linear(numcelltran,numcell,numendcell,maxpct,dxc)
%% Maarten Buijsman, GFDL, 2011-09-29
%% gridstretching, with gradual increase to maxpct
%%
%% numcell=number of cells in stretched grid
%% numcelltran=transition from 0 to maxpct
%% numendcell=number of cells + 1 that have a constant value
%% maxpct=max percentage
%% dxc=undisturbed center grid cell

function dxd = stretch_grid_linear(numcelltran,numcell,numendcell,maxpct,dxc);

if numcelltran>numcell; disp('ERROR'); return; end

%% test
% dxc = 2000;
% numcell     = 48;
% numcelltran = 10;
% numendcell  = 12;
% maxpct      = 0.036;


dxd=ones(1,numcell);
dxs = dxc;
for i=1:numcell-numendcell
    dxd(i) = dxs + dxs*min(maxpct,maxpct*i/numcelltran);
    dxs = dxd(i);
    
    maxpcti(i) = min(maxpct,maxpct*i/numcelltran);
end
dxd(numcell-numendcell:end) = dxs;

return

%% test
% figure
% plot(dxd,'k.-')