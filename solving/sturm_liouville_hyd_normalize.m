%% function [C,Cg,L,Weig2,Ueig2] = sturm_liouville_hyd_normalize(om,Nb,DZ,f)
%% Maarten Buijsman, GFDL, 2012-05-29
%% solves non-hydrostatic sturm liouville equation
%% input stratification N(z) (bottom up) at cell faces, 
%% DZ, Coriolis frequency f, and omega om (rad/s)
%% computes normalized eigen functions Weig and Ueig (=Peig) at cell centers
%% 0th mode is barotropic mode
%% computes phase speed C, group speed Cg, and wave length L

function [C,Cg,L,Weig2,Ueig2] = sturm_liouville_hyd_normalize(om,Nb,DZ,f)

% %% test =========================================
% clear all
% dirout = '/home/m2b/data/projects/SCS_batan/';
% fnm    = 'rho_merged_alford.mat';
% load([dirout fnm])
% H = 1000;
% 
% f1=figure
% subplot(1,2,1)
% plot(rhoi2,zi2)
% holder
% 
% %% extrapolate linearly
% dz1 = 10;
% zw   = [-H:dz1:0];  %bottom up !!!!!!
% 
% %% compute N
% gravity=9.81; rhoNil=999.8;
% N2r = -gravity/rhoNil*diff(rhoi2)./diff(zi2);             % at rho
% zi3 = zi2(1:end-1)/2 + zi2(2:end)/2;
% N2  = [interp1(zi3,N2r,zw(1:end-1),'linear','extrap') 0]; % at zw
% 
% Nbr = sqrt(N2);
% %Nbr = sqrt(N2)*0+1e-2;
% 
% figure(f1)
% subplot(1,2,2)
% plot(Nbr,zw)
% 
% om = 12.1408331767/24/3600;
% z2 = zw; 
% Nb = Nbr;
% f  = 1e-5;
% 
% DZ = dz1;
% %% test =========================================

%% make sure Nb is column vector
[a,b]=size(Nb);
if a<b; Nb=Nb'; end

N  = length(Nb)-1; %% number of layers
H  = N*DZ;

%% generalized + alternative
A=[]; B=[];
i=1;
A(i,i)   = -2*1/(DZ)^2;
A(i,i+1) =  1*1/(DZ)^2;
for i=2:N-2
    A(i,-1+i) =  1*1/(DZ)^2;
    A(i, 0+i) = -2*1/(DZ)^2;
    A(i, 1+i) =  1*1/(DZ)^2;    
end
i=i+1;
A(i,-1+i) =  1*1/(DZ)^2;
A(i, 0+i) = -2*1/(DZ)^2;

%% hydrostatic part
N2 = Nb.^2; 
B = diag(-N2(2:end-1))/(om^2-f^2);  % N to power 2

%% solve it ===========================================
[W1,k2] = eig(A,B);

%% wave lengths and phase speeds
k = sqrt(diag(k2));
[k,Is]=sort(k,'ascend'); %sort
W1 = W1(:,Is);
C  = om./k;
L  = 2*pi./k;
Cg = (om^2-f^2)./(om*k);


%% get U/P structure functions
W2    = [zeros(1,size(W1,1)); W1; zeros(1,size(W1,1))];
Ueig1 = [ones(N,1) (W2(2:end,:)- W2(1:end-1,:))/DZ];


%% Normalize Ueig2
%AA=repmat(nansum(Ueig1.^2.*DZ,1)./H,[N 1]).^(1/2);
AA=repmat(sum(Ueig1.^2.*DZ,1)./H,[N 1]).^(1/2);
AA(AA==0)=Inf;
Ueig2=Ueig1./AA;
%Ueig2(:,Ueig2(N,:)<0)=-Ueig2(:,Ueig2(N,:)<0);
Ueig2(:,Ueig2(1,:)<0)=-Ueig2(:,Ueig2(1,:)<0); %% bottom

%% get W structure functions
Weig1 = [zeros(N,1) W2(1:end-1,:)/2 + W2(2:end,:)/2]; %at cell centers

% SS = size(W2(:,1),1);
% figure; 
% plot(W2(:,1),1:SS,'k.-')
% hold
% plot(Weig1(:,2),[1:SS-1]+0.5,'r.-')

%% normalize Weig2 (why????)
N2c = N2(1:end-1)/2 + N2(2:end)/2; % at cell centers
AA=repmat(sum(Weig1.^2.*repmat(N2c,[1 N])*DZ,1)./H,[N 1]).^(1/2);
%AA=repmat(sum(Weig1.^2.*DZ,1)./H,[N 1]).^(1/2);
%bad? AA=repmat(sum(Ueig1.^2.*DZ,1)./H,[N 1]).^(1/2);
AA(AA==0)=Inf;
Weig2=Weig1./AA;
%Weig2(:,Weig2(N,:)<0)=-Weig2(:,Weig2(N,:)<0);
Weig2(:,Weig2(1,:)<0)=-Weig2(:,Weig2(1,:)<0); %% bottom

%disp(AA(1,:))

% below is Kelly
% A=repmat(nansum(PHI2.^2.*repmat(N2,[1 N])*dz,1)./H,[N 1]).^(1/2);
% A(A==0)=Inf;
% PHI2=PHI2./A;
% PHI2(:,PHI2(N,:)<0)=-PHI2(:,PHI2(N,:)<0);


%figure; plot(Ueig2(:,1:10),1:N)
%figure; plot(Weig2(:,1:10),1:N)
%sum(Ueig2(:,2).*Ueig2(:,3),1)
%sum(Ueig2(:,10).*Ueig2(:,10)*DZ,1)/H

