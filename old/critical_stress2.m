%% critical_stress.m
%% Maarten Buijsman, NIOZ, 25-10-05
%% calculates critical shear stress and shiels param.

clear all

grav = 9.81;
rho_s = 2650;
rho = 1027;
%d50 = 300/1e6;
d50 = [200:100:600]/1e6;
%d50 = [1:1:10000]/1e6;
Temp = 10;

%Temp = [0:1:20]; 
%visc = 1/10*(10.^(1301./(998.333+8.1855.*(Temp-20)+0.00585.*(Temp-20).^2)-3.30233))/1000; %% fresh
visc = 1/10*(10.^(1301./(998.333+8.1855.*(Temp-20)+0.00585.*(Temp-20).^2)-3.30233))/1000+0.05*1e-6; %% salt
%plot(Temp,visc)
%visc_i = interp1(Temp,visc,10);

%% Soulsby, 97, eq77
s = rho_s/rho;
Dstar = (grav*(s-1)/visc^2)^(1/3)*d50;
theta_cr = 0.30./(1+1.2*Dstar)+0.055.*(1-exp(-0.020*Dstar));
tau_cr = theta_cr.*grav*(rho_s-rho).*d50;
Dstar = (grav*(s-1)/visc^2)^(1/3)*d50;

%% tau
figure
plot(d50*1e6,tau_cr,'r')
%loglog(d50*1e6,tau_cr,'r')
hold;
plot_label02(10,0.5,'\tau_{cr}(T)','d50 [mu]','\tau_{cr} [Nm^{-2}]')
legend('0','10','20');

%% theta
figure
Dstar
loglog(Dstar,theta_cr,'r')
%plot(d50*1e6,theta_cr,'r')
% hold;
% plot_label02(10,0.5,'\theta_{cr}(T)','d50 [mu]','\theta_{cr} [-]')
legend('0','10','20');

%%-----------------------------------------------------------------------
%% shear stress due to currents
%%-----------------------------------------------------------------------
uc = 1;
depth=17;
tau = rho*grav*(uc./(18*log10(12*depth./(2.5*d50)))).^2;

% ustar = 9.02/100;
% tau = ustar^2*1027

%% soulsby p54
%% various ways to calculate skin friction 
%% ann all about the same order....
C100 = 0.0024;
d50 = 400/1e6;
Uave = 1.25;
h=25;
U100=0.8;
alpha=0.0475; beta=2/7;

%% 1
ustar = 1/7*(d50/h)^(1/7)*Uave;
tau_c = rho*ustar^2;

%% 2
tau_c = rho*C100*U100^2;

%% 3
zo = d50/12;
Cd = alpha*(zo/h)^beta;
%% or
Cd = (0.40/(1+log(zo/h)))^2;
tau_c = rho*Cd*Uave^2;

%% van Rijn's method
tau_c = rho*grav*(Uave/(18*log10(12*h/(2.5*d50))))^2;
%%-----------------------------------------------------------------------
%% due to waves
%%-----------------------------------------------------------------------
%h=25;
h=20;
Tw = 4;
Hs = 2;

% h=10;
% Tw = 8;
% Hs = 3;


om = 2*pi/Tw;

%% estimate wave number k
%% Soulsby 97, p71
epsi = om^2*h/grav;
if epsi>1
    eta = epsi*(1+0.2*exp(2-2*epsi));
else
    eta = epsi^(1/2)*(1+0.2*epsi);
end
k=eta/h;
L = 2*pi/k;

%% get orbital velocity near bottom
Hw=Hs/sqrt(2);
%Hw = 0.05;

Uw = pi*Hw/(Tw*sinh(k*h));

%%-----------------------------------------------------------------------
%% due to waves and currents
%%-----------------------------------------------------------------------

clear all

grav = 9.81;
rho_s = 2650;
rho = 1027;

%% table 10, p95
%% fw, friction factor
Azo = [1e2 1e3 1e4 1e5];
GM79Azo = [0.1057 0.0316 0.0135 0.00690];
Azo_x = [1e1 Azo 1e6 1e7];

figure
plot(log10(Azo),log10(GM79Azo))

%% 1
[x_s,y_fit,r,p,cf_fw1] = line_fit(log10(Azo),log10(GM79Azo),1);
hold; plot(x_s,y_fit,'r--')

X = log10(Azo_x);
Y1Azo = polyval(cf_fw1,X);
plot(X,Y1Azo,'k--')

%% 2
[x_s,y_fit,r,p,cf_fw2] = line_fit(log10(Azo),log10(GM79Azo),2);

X = log10(Azo_x);
Y2Azo =  polyval(cf_fw2,X);
plot(X,Y2Azo,'c--')

figure;
plot(Azo,GM79Azo,'r--')
hold
plot(Azo_x,10.^Y1Azo,'k--')
plot(Azo_x,10.^Y2Azo,'c--')

%% Cd, drag coefficient
zoh = [1e-2 1e-3 1e-4 1e-5];
GM79zoh = [0.01231 0.00458 0.00237 0.00145];
zoh_x = [1e-1 zoh 1e-6 1e-7];

figure
plot(log10(zoh),log10(GM79zoh))

%% 1
[x_s,y_fit,r,p,cf_Cd1] = line_fit(log10(zoh),log10(GM79zoh),1);

%% 2
[x_s,y_fit,r,p,cf_Cd2] = line_fit(log10(zoh),log10(GM79zoh),2);
hold; plot(x_s,y_fit,'r--')

X = log10(zoh_x);
Y2zoh = polyval(cf_Cd2,X);
plot(X,Y2zoh,'c--')

figure;
plot(zoh,GM79zoh,'r--')
hold
plot(zoh_x,10.^Y2zoh,'c--')

%% BEGIN HERE
% Soulsby example
% Uw = 0.5;
% T=12.6;
% phi = 30*pi/180;
% Uave = 1;
% zo = 0.001;
% h=10;
% rho=1027;

Uw = 0.04;
T=4;
phi = 0*pi/180;
Uave = 1.5;
d50 = 400/1e6;
%zo = 0.001;
zo = d50/12;
h=25;
rho=1027;

%% functions table 10
%% calc the mean of linear and quadratic fit
%% fw 1
A = Uw*T/(2*pi);
x = log10(A/zo);
fw = 0; fc = cf_fw1;
for i=1:length(fc)
    fw = fw + fc(i)*x^(length(fc)-i);
end
fw1 = fw;

%% fw 2
fw = 0; fc = cf_fw2;
for i=1:length(fc)
    fw = fw + fc(i)*x^(length(fc)-i);
end
fw2 = fw;
fw = 10^((fw1+fw2)/2);

%% Cd 1
x = log10(zo/h);
Cd = 0; fc = cf_Cd1;
for i=1:length(fc)
    Cd = Cd + fc(i)*x^(length(fc)-i);
end
Cd1 = Cd;

%% Cd 2
Cd = 0; fc = cf_Cd2;
for i=1:length(fc)
    Cd = Cd + fc(i)*x^(length(fc)-i);
end
Cd2 = Cd;
Cd = 10^((Cd1+Cd2)/2); %% about the same as from equations above...

%% current only
tau_c = rho*Cd*Uave^2;

%% wave only
tau_w = 1/2*rho*fw*Uw^2;

X = tau_c/(tau_c+tau_w);

%% coefficients
%% table 9, p95
as = [0.11 1.95 -0.49 -0.28];
ms = [0.65 -0.22 0.15 0.06];
ns = [0.71 -0.19 0.17 -0.15];
I = 0.67;
bs = [0.73 0.40 -0.23 -0.24];
ps = [-0.68 0.13 0.24 -0.07];
qs = [1.04 -0.56 0.34 -0.27];
J = [0.5];

%% and formulas...
a = (as(1)+as(2)*abs(cos(phi))^I) + (as(3)+as(4)*abs(cos(phi))^I)*log10(fw/Cd);
m = (ms(1)+ms(2)*abs(cos(phi))^I) + (ms(3)+ms(4)*abs(cos(phi))^I)*log10(fw/Cd);
n = (ns(1)+ns(2)*abs(cos(phi))^I) + (ns(3)+ns(4)*abs(cos(phi))^I)*log10(fw/Cd);

b = (bs(1)+bs(2)*abs(cos(phi))^J) + (bs(3)+bs(4)*abs(cos(phi))^J)*log10(fw/Cd);
p = (ps(1)+ps(2)*abs(cos(phi))^J) + (ps(3)+ps(4)*abs(cos(phi))^J)*log10(fw/Cd);
q = (qs(1)+qs(2)*abs(cos(phi))^J) + (qs(3)+qs(4)*abs(cos(phi))^J)*log10(fw/Cd);

Z =    1+a*X^m*(1-X)^n;
Y = X*(1+b*X^p*(1-X)^q);

tau_max = Z*(tau_c+tau_w);
tau_m   = Y*(tau_c+tau_w);




