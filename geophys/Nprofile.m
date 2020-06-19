%% function Nw = Nprofile(zr,rhor,zw,grav,rho0,plotfig)
%% Maarten Buijsman, UCLA, 
%% compute stratification
%% vectors zr rhor at rho-points
%% vectors Nw and zw at w-pints, zw should include bottom!

function Nw = Nprofile(zr,rhor,zw,grav,rho0,plotfig);

% plotfig = true
% zr   = zz2(:,1);
% rhor = rho_rest(:,1);
% zw = zz2w(:,1);

if min(zw) > 0; disp(['z should be positive!']); return; end

zi   = linspace(min(zw),0,501);
%rhoi = interp1(zr,rhor,zi,'cubic'); % extrapolate to surface/bottom
rhoi = interp1(zr,rhor,zi,'linear','extrap'); % extrapolate to surface/bottom

zi2 = (zi(1:end-1) + zi(2:end))/2;
N2  = sqrt(-grav/rho0*diff(rhoi)./diff(zi));
%Nw  = interp1(zi2,N2,zw,'cubic');  % extrapolate to surface/bottom
Nw  = interp1(zi2,N2,zw,'linear','extrap');  % extrapolate to surface/bottom

if plotfig==true
    figure;
    subplot(2,1,1);plot(rhor,zr,'r-',rhoi,zi,'b--')
    subplot(2,1,2);plot(N2,zi2,'r-',Nw,zw,'b--')    
end
