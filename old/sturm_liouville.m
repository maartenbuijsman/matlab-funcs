%% function [Cf,Lf,C,L,gam,V] = sturm_liouville(om,z2,Nb,f)
%% computes phase speed C and wave length L
%% and eigen values and vectors gam and V
%% given stratification N(z), Coriolis forcing f, and omega om (rad/s)
%% frequency independent phase speed is C = 1./sqrt(-gam);
%% Cf and Lf are with Coriolis forcing
%% all values have to be bottom up

function [Cf,Lf,C,L,gam,V] = sturm_liouville(om,z2,Nb,f)

%% test
%z2=zw; Nb=NS;
flg= 0;
le = length(z2);
DZ = mean(diff(z2));
botflg = 0;

if DZ<0;               
    %% reverse
    z2 = z2(end:-1:1);
    Nb = Nb(end:-1:1);
    DZ = mean(diff(z2));
    botflg=1;   
end
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

for i=1:length(z2)-2
    B(i,i) = Nb(i+1)^2; % N to power 2
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

%% reverse again
if botflg
    V = V(end:-1:1,:);
end

% same result
%C=sqrt(-1./gam);
%L=NaN;

