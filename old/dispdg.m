function [y,eta,k1,omega]=dispdg(h,d,f,C,eigv);
%h=pos. depth prof
%d=dx
%f=cor freq
%C=phase speed
%eigv = eigen values

% % h=fliplr(Di);
% h=Di;
% d=d;
% f =1.0423e-04;
% C=200;
% eigv=100;

clear k2;
H=max(h);
delxm=d;
delxkm=d/1000;
g = 9.8;
N = length(h);
% spectral parameter 

% spectral parameter 
%C = 30; % phase speed = omega/k [m/s]
S = f/C;
% coefficients A and B are specified:
for i=1:N
    A(i) = g*h(i)*S*S/(S*S*g*h(i)-f*f);
end

for i=2:N-1
    B(i)=(h(i+1)-h(i-1))/(2*delxm*h(i));
end
B(1)=(h(2)-h(1))/(delxm*h(1));
B(N)=(h(N)-h(N-1))/(delxm*h(N));

% Define matrix coefficients
aa(1:N,1:N) = 0;
for i=2:N-1
    aa(i,i-1) = A(i)*(1/(delxm*delxm)-B(i)/(2*delxm)); 
    aa(i,i) = A(i)*(S*B(i)-f*f/(g*h(i))-2/(delxm*delxm));
    aa(i,i+1) = A(i)*(1/(delxm*delxm)+B(i)/(2*delxm));
end

%1st line:
aa(1,1) = 2*A(1)*S/delxm-2*A(1)/(delxm*delxm)-A(1)*f*f/(g*h(1));
aa(1,2) = 2*A(1)/(delxm*delxm);

%last (Nth) line:
aa(N,N-1) = 2*(A(N)*S - A(N)*B(N))/(S*delxm*delxm+2*delxm);
aa(N,N) = A(N)*B(N)*S-A(N)*f*f/(g*h(N))+2*(A(N)*B(N)-A(N)*S)/(S*delxm*delxm+2*delxm); 

%calculate eigenvalues
ksq=eig(aa);

%calculate particular wave structure (eigenvector)
% k2 = ksq(161);
%find 0 mode
% mod0=find(ksq>min(abs(ksq)),1,'last');
display(find(ksq>min(abs(ksq)),3,'last'));
%calculate particular wave structure (eigenvector)
kk = ksq(eigv);
% if kk>0
    k2=kk;
% wave properties below are nor used in calculations, just diagnostics
k1=sqrt(k2);
lmbd=2*pi/(k1*1000);
omega=C*k1;
%omega=-C*k1; % we use it for negative phase speed
period=2*pi/(omega*3600);
%------------------------------------------------------------------------
eta(1:N) = 1;
eta(2) = (k2-aa(1,1))*eta(1)/aa(1,2);
for i=2:N-1
    eta(i+1) = (eta(i)*(k2-aa(i,i))-eta(i-1)*aa(i,i-1))/aa(i,i+1);
end

for i=1:N
    y(i)=delxkm*(i-1);
end
figure (2)
subplot(2,1,1);
plot(y,eta);
eta1=eta';
%save mode0_125.dat eta1 -ascii;
xlabel('offshore distance [km]');
ylabel('sea level elevation [m]');
%axis([0 300 -2 2]);
%title('C=120 m/s, zero mode, wavelength=5081 km, T=11.76 hrs');
end