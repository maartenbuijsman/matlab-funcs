%% wave_shear_stress
%% 

clear all

grav=9.81;
h=20;
Tz = 4; tab = 0.05
%Tz = 5; tab = 0.09
Hs = 2;

rho=1025;
d90 = 500*1e-6 % meter

% h=10;
% Tz = 8;
% Hs = 3;

om = 2*pi/Tz;

%% estimate wave number k
%% Soulsby 97, p71
epsi = om^2*h/grav;
if epsi>1
    eta = epsi*(1+0.2*exp(2-2*epsi));
else
    eta = epsi^(1/2)*(1+0.2*epsi);
end
k =eta/h;
L = 2*pi/k;

%% check if shallow or deep
hch = 0.1*grav*Tz^2;
if h<hch; disp(['shallow; ',num2str([h hch])]); 
else disp('not shallow and little efect waves'); end

%% get orbital velocity near bottom
Hw=Hs/sqrt(2);
Uw = pi*Hw/(Tz*sinh(k*h));

%% from figure 14 Soulsby
%h=10; Hs=3; Tz=8;
Tn=sqrt(h/grav);
Tn/Tz
%tab = 0.205
Urms = tab*Hs/Tn


%% estimation of skin tau_w 
Uw = sqrt(2)*Urms;
T  = Tz*1.281;

ks=3*d90;
z0 = ks/30;
%z0=6*1e-3;
A = Uw*T/(2*pi);
fw = 1.39*(A/z0)^(-0.52);
tau_w = 0.5*rho*fw*Uw^2;

%% compare with current
uc=1.25; 

tau_c = rho*grav*(uc./(18*log10(12*h./(3*d90)))).^2;







