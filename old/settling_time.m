%% settling_time.m
%% Maarten Buijsman, 14-09-07
%% play with fall velocity and settling time

fact = 0.66;
d50 = 350*1e-6;
dp = d50*fact;
grav = 9.81; 
rho = 1025;
visc = 1e-6;
kappa=0.4;
h2=25;
h=25;

s = 2650./rho; 
Dstar = ((s-1).*grav./visc.^2).^(1/3)*dp;
ws = visc/dp.*(sqrt(10.36^2+1.049*Dstar.^3)-10.36);

%ws = 0.023 %in winter

T = h2/ws/60 %settling time

%u=0.4;
%D = T*u*60

ks=1.1*0.5*(1-exp(-25*5/7));
C = 18*log10(12*h/ks);

%ustar = sqrt(grav)*U/C;
Z = 1;
%beta = 1+2*(ws/ustar)^2;
beta=2;
ustar=ws/(beta*kappa*Z);

%% velocity when ustar=ws;
%ustar=ws;
U=ustar/sqrt(grav)*C;

D = T*U*60 %horizontal distance
