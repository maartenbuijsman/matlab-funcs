%% function [Cf,Lf,C,L,gam,V] = sturm_liouville_nonhyd_Ndr(om,z2,Nb,f)
%% Maarten Buijsman, UCLA, 2009-05-04
%% solves non-hydrostatic sturm liouville equation
%% computes phase speed C and wave length L
%% and eigen values and vectors gam and V
%% given stratification N(z), Coriolis forcing f, and omega om (rad/s)
%% Cf and Lf are with Coriolis forcing
%% all values have to be bottom up
%% sturm_liouville_nonhyd_Ndr.m => as Kelly et al (2012) d/dz[(1/(N2-om2) dW/dz] + (1-f2/om2)/cg2 W

%% something fishy here .......
%% vertical velocity does not correctly translate in horiz velocity when vertical gradient is taken .....

%% HA - in Kelly 2012 the eigenfunction problem solves for the horizontal velocity and pressure ....
%% sooooo, this is not the vertical vlocity ...... 


function [Cf,Lf,C,L,gam,V] = sturm_liouville_nonhyd_Ndr(om,z2,Nb,f)

% % %% test =========================================
% % clear all
% % dirout = '/home/m2b/data/projects/SCS_batan/';
% % fnm    = 'rho_merged_alford.mat';
% % load([dirout fnm])
% % H = 1000;
% % 
% % f1=figure
% % subplot(1,2,1)
% % plot(rhoi2,zi2)
% % holder
% % 
% % %% extrapolate linearly
% % dz1 = 10;
% % zw   = [-H:dz1:0];  %bottom up !!!!!!
% % 
% % %% compute N
% % gravity=9.81; rhoNil=999.8;
% % N2r = -gravity/rhoNil*diff(rhoi2)./diff(zi2);             % at rho
% % zi3 = zi2(1:end-1)/2 + zi2(2:end)/2;
% % N2  = [interp1(zi3,N2r,zw(1:end-1),'linear','extrap') 0]; % at zw
% % 
% % Nbr = sqrt(N2);
% % %Nbr = sqrt(N2)*0+1e-2;
% % 
% % figure(f1)
% % subplot(1,2,2)
% % plot(Nbr,zw)
% % 
% % om = 12.1408331767/24/3600;
% % z2 = zw; 
% % Nb = Nbr;
% % f  = 1e-5;
% % %% test =========================================

% om = omfq;
% z2=zw; 
% Nb=Nbfq; 
% f=fcor;

flg= 0;
le = length(z2);
DZ = mean(diff(z2));

if DZ<0;                disp(['ERROR, z and N need to be bottom up!!!']); return; end
if std(diff(z2))>0.001; disp(['   INTERPOLATE, z is not equidistant!!!']); 
   %% interpolate on constant grid
   ztmp=z2; z2=[];    
   Ntmp=Nb; Nb=[];
   z2 = linspace(ztmp(1),ztmp(end),100);
   Nb = interp1(ztmp,Ntmp,z2,'cubic'); 
   DZ = mean(diff(z2)); 
   flg=1;
end

% % %% OLD ==================================================================
% % %% generalized + alternative
% % A=[]; B=[];
% % i=1;
% % A(i,i)   = -2*1/(DZ)^2;
% % A(i,i+1) =  1*1/(DZ)^2;
% % for i=2:length(z2)-3
% %     A(i,-1+i) =  1*1/(DZ)^2;
% %     A(i, 0+i) = -2*1/(DZ)^2;
% %     A(i, 1+i) =  1*1/(DZ)^2;    
% % end
% % i=i+1;
% % A(i,-1+i) =  1*1/(DZ)^2;
% % A(i, 0+i) = -2*1/(DZ)^2;
% % 
% % %% non-hydrostatic part: Nb-omega 
% % N2 = Nb.^2-om^2; Isel=find(N2<0); 
% % 
% % if length(Isel)>0; 
% % %    disp(['omega2 ' num2str(om) ' > Nb; number = ' num2str(length(Isel))]); 
% % %    figure; plot(Nb,z2,'k-',Nb(Isel),z2(Isel),'rs')
% %     Nb(Isel) = om+0.0001e-4; %% to avoid singularities
% % end
% % 
% % for i=1:length(z2)-2
% %     B(i,i) = Nb(i+1)^2-om^2;  % N to power 2
% % end
% % %% OLD ==================================================================

%% non-hydrostatic part: Nb-omega 
Nb(end) = Nb(end-1)/2+Nb(end)/2;
N2 = Nb.^2-om^2; Isel=find(N2<0); 

%% how do we deal with negative or zero N2 values .....
if length(Isel)>0; 
%    disp(['omega2 ' num2str(om) ' > Nb; number = ' num2str(length(Isel))]); 
%    figure; plot(Nb,z2,'k-',Nb(Isel),z2(Isel),'rs')
    N2(Isel) = 1e-12; %% to avoid singularities
    
    if ~isempty(Isel)
        disp(['singularities at depth .... ' num2str(Isel) ' of ' num2str(length(z2))])
    end
end

N2inv = 1./N2;

%% generalized + alternative
%% comments page 225-228, 2012/05/25
A=zeros(length(z2)-2,length(z2)-2); B=A;

%% at bottom
i=1;
A(i,i)   =  - N2inv(i) - 2*N2inv(i+1) - N2inv(i+2);
A(i,i+1) =                 N2inv(i+1) + N2inv(i+2);
for i=2:length(z2)-3  %% z=3!!!!!!!!!!!!!!
    A(i,-1+i) =   N2inv(i) +   N2inv(i+1);
    A(i, 0+i) = - N2inv(i) - 2*N2inv(i+1) - N2inv(i+2);
    A(i, 1+i) =                N2inv(i+1) + N2inv(i+2);    
end

%% at surface
i=i+1;
A(i,-1+i) =   N2inv(i)  +   N2inv(i+1);
A(i, 0+i) = - N2inv(i)  - 2*N2inv(i+1) - N2inv(i+2);

A = A/(2*DZ^2);

for i=1:length(z2)-2
    B(i,i) = 1;  % N to power 2
end

%% solve it ===========================================
[V,gam2] = eig(A,B);
%V2=[]; V2(2:size(V,1)+1,:) = V; V2(end+1,:) = 0;
gam = diag(gam2);

% with coriolis force
k = sqrt(-gam*(om^2-f^2));
Cf =om./k;
Lf = 2*pi./k;

% without coriolis force
f=0;
k = sqrt(-gam*(om^2-f^2));
C =om./k;
L = 2*pi./k;

%% map to org. grid
if flg==1;
    Vtmp = V;
    V = interp1(z2(2:end-1)',Vtmp,ztmp(2:end-1),'cubic');
%    figure; plot(Vtmp(:,end-1),z2(2:end-1)','k-',V(:,end-1),ztmp(2:end-1),'b-')
end

% % %% test
% % figure
% % subplot(3,1,1); 
% % plot(V(:,end),z2(2:end-1),'b.-')
% % hold
% % plot(V(:,end-1),z2(2:end-1),'r.-')
% % plot(V(:,end-2),z2(2:end-1),'g.-')
% % title('eig')
% % 
% % subplot(3,1,2); 
% % plot(1:length(z2)-2,Cf,'b.-')
% % title('Cf')
% % 
% % subplot(3,1,3); 
% % plot(1:length(z2)-2,Lf/1e3,'b.-')
% % title('Lf')


% same result
%C=sqrt(-1./gam);
%L=NaN;

