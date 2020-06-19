%% un_dro.m
%% balance for paper

rho1 = 0.75;rho2 = 2.25;
%rho1 = 0.75;rho2 = 10.25;
grav = 9.81;
rho = (rho1+rho2)/2;
B = 4000; H=20;
W=B/2; h=H/2;
T = 12.5*3600;
DT = [0:100:T/2];
un = 0.05;

DPDN = grav/rho*(rho2-rho1)*un*h/W^2*DT;

plot(DT,DPDN*1e4,'k.')

Dh = 300*un*h/W