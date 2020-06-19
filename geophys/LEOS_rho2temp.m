%% function T=LEOS_rho2temp(rho,rho0,T0,alp)
%% MCB, UCLA, 2009-03-02
%% lin. eq. of state
%% input pot. density, output pot. temp
%% alp is therm exp. coef. 

function T=LEOS_rho2temp(rho,rho0,T0,alp)

T = T0+1/alp*(1-rho/rho0);

return

% %% rho is computed like this, here
% T - T0 = +1/alp*(1-rho/rho0);
% (T - T0)*alp = 1-rho/rho0;
% (T - T0)*alp - 1 = -rho/rho0;
% -((T - T0)*alp - 1)*rho0 = rho;
% 
% rho = rho0-((T - T0)*alp*rho0
% 
% %% rho is computed like this, in ROMS
% rho=rho0 -alp*(T)