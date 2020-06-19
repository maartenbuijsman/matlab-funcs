%% function rho=LEOS_temp2rho(T,rho0,T0,alp)
%% MCB, UCLA, 2009-03-02
%% lin. eq. of state
%% input pot. T, output pot. density
%% alp is therm exp. coef. 

function rho=LEOS_temp2rho(T,rho0,T0,alp)

rho = rho0*(1-alp*(T-T0));      

return
