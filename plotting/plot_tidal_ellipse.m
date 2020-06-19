%% function [u,v] = plot_tidal_ellipse(au,bu,av,bv,fq,time_d);
%% MCB, UCLA, 2008-03-28
%% plot tidal ellipse, includes line where t=0;
%% note that, when comparing phases (from t=0 to max u) by plotting ellipses
%% the direction of rotation is important (imagine reversing ellipse so rot. direction is the same) !!!!

function [u,v] = plot_tidal_ellipse(au,bu,av,bv,fq);

%% tidal ellipse
ts = 0; 
t = linspace(ts,ts+2*pi/fq,50);    %time in days
u1 = au*cos(fq*t) + bu*sin(fq*t); 
v1 = av*cos(fq*t) + bv*sin(fq*t); 

%% this is not phase
%% include phase
u = [0 u1];
v = [0 v1];


return



%% test 1 -------------------------------
a = 1;
time_d =[0:25]/24;
fq = 12.1408331767;
ud =a*cos(fq*time_d-45*pi/180);
vd =a/3*cos(fq*time_d-30*pi/180);

plot_fg=0;
freq_sel=[46];
[Ao_ud,au,bu,fq2,coef_det,sigma,data,fit] = harmonic03(time_d,ud,freq_sel,1,plot_fg,freq_dat); 
[Ao_vd,av,bv,fq2,coef_det,sigma,data,fit] = harmonic03(time_d,vd,freq_sel,1,plot_fg,freq_dat); 
%atan(bv/av)*180/pi

[Amax,Amin,Ec,phi,psi,Rplus,Rmin,phiplus,phimin]=tidal_ellipse2(au,bu,av,bv);    
phase = phi*180/pi
inclin= psi*180/pi
phase-inclin
tmax = phi/fq

%% tidal ellipse
ts = 0; %time in days
%ts = time_d(1); %time in days
t = linspace(ts,ts+2*pi/fq,50);
u1 = au*cos(fq*t) + bu*sin(fq*t); 
v1 = av*cos(fq*t) + bv*sin(fq*t); 

%% this is not phase
%% include phase
u = [0 u1];
v = [0 v1];

f1 = figure;
s1=subplot(2,1,1);
plot(u,v,'k.-')
for i=1:length(u); text(u(i),v(i),num2str(i-1),'verticalalignment','middle','fontsize',6); end
%hold on; plot([-1.1 1.1],tan(psi)*[-1.1 1.1])
axis equal

s2=subplot(2,1,2);
plot(t,sqrt(u(2:end).^2+v(2:end).^2),'k.-')

t(7)


%% test 2 -------------------------------
a = 1.2;
ud =a*cos(fq*time_d-50*pi/180);
vd =a/3*cos(fq*time_d-20*pi/180);

[Ao_ud,au,bu,fq2,coef_det,sigma,data,fit] = harmonic03(time_d,ud,freq_sel,1,plot_fg,freq_dat); 
[Ao_vd,av,bv,fq2,coef_det,sigma,data,fit] = harmonic03(time_d,vd,freq_sel,1,plot_fg,freq_dat); 
%atan(bv/av)*180/pi

[Amax,Amin,Ec,phi,psi,Rplus,Rmin,phiplus,phimin]=tidal_ellipse2(au,bu,av,bv);    
phase2 = phi*180/pi
inclin2= psi*180/pi
phase2-inclin2

tmax = phi/fq
%% test -------------------------------


%% tidal ellipse
u1 = au*cos(fq*t) + bu*sin(fq*t); 
v1 = av*cos(fq*t) + bv*sin(fq*t); 

%% this is not phase
%% include phase
u = [0 u1];
v = [0 v1];

%% test -------------------------------
figure(f1);
subplot(s1);
hold on; plot(u,v,'r.-')
for i=1:length(u); text(u(i),v(i),num2str(i-1),'verticalalignment','middle','fontsize',6,'color','r'); end
axis equal
hold off;

subplot(s2);
hold on; plot(t,sqrt(u(2:end).^2+v(2:end).^2),'r.-')
hold off;

phase2-phase














