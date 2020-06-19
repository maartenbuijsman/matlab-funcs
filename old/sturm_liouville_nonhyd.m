%% function [Cf,Lf,C,L,gam,V] = sturm_liouville_nonhyd(om,z2,Nb,f)
%% Maarten Buijsman, UCLA, 2009-05-04
%% solves non-hydrostatic sturm liouville equation
%% computes phase speed C and wave length L
%% and eigen values and vectors gam and V
%% given stratification N(z), Coriolis forcing f, and omega om (rad/s)
%% Cf and Lf are with Coriolis forcing
%% all values have to be bottom up

function [Cf,Lf,C,L,gam,V] = sturm_liouville_nonhyd(om,z2,Nb,f)

%% test
om = 2*pi/(12*3600);
f=0;
z2 = [-1000:10:0];
Nb = ones(size(z2))*0.001;

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

%% generalized + alternative
A=[];
i=1;
A(i,i)   = -2*1/(DZ)^2;
A(i,i+1) =  1*1/(DZ)^2;
for i=2:length(z2)-3
    A(i,-1+i) =  1*1/(DZ)^2;
    A(i, 0+i) = -2*1/(DZ)^2;
    A(i, 1+i) =  1*1/(DZ)^2;    
end
i=i+1;
A(i,-1+i) =  1*1/(DZ)^2;
A(i, 0+i) = -2*1/(DZ)^2;

%% non-hydrostatic part: Nb-omega 
N2 = Nb.^2-om^2; Isel=find(N2<0); 

if length(Isel)>0; 
%    disp(['omega2 ' num2str(om) ' > Nb; number = ' num2str(length(Isel))]); 
%    figure; plot(Nb,z2,'k-',Nb(Isel),z2(Isel),'rs')
    Nb(Isel) = om+0.0001e-4; %% to avoid singularities
end

for i=1:length(z2)-2
    B(i,i) = Nb(i+1)^2-om^2;  % N to power 2
end
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

% same result
%C=sqrt(-1./gam);
%L=NaN;

%figure; plot(V(:,end),z2(2:end-1)','k-')
