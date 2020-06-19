%% plot_tidebeam
%% MCB, UCLA, 2008-09-23
%% gives output of beams at one location
%% Isel has the indices of x-coordinate

function [xpos,ypos]=plot_tidebeam(omega,fnam1,Isel);

%% get directory
dir1 = '/work/WORK/experiments/'

%% load data
D = dir(dir1); %% search directory for contents

%% macth patterns and fill fnames
j=0; j2=0;
for i=1:length(D)
    x = strmatch([fnam1],D(i).name); if x==1; j=j+1;   fnames(j)   = {D(i).name}; end
end

%% get h values
nc=netcdf([dir1 char(fnames(1))]);
dt = getfield(nc(2),'ocean_time')/(getfield(nc(2),'time_step')-1); %get time step
u_u=squeeze(nc{'u'}); hi=size(u_u,2); wi=size(u_u,3); le=size(u_u,4); clear u_u

%% get z and x values
Itime=1; Iy=2;
zeta=squeeze(nc{'zeta'}(Itime,:,:)); h=squeeze(nc{'h'}(:,:));
theta_s=nc.theta_s(:); theta_b=nc.theta_b(:); hc=nc.hc(:);
N=hi;type='r'; zz=zlevs(h,zeta,theta_s,theta_b,hc,N,type); zz2 = squeeze(zz(:,Iy,:));
xgrid=squeeze(nc{'x_rho'}(:,:)); xx = ones(size(zz2,1),1)*xgrid(1,:); xg = xgrid(1,:);
pm = squeeze(nc{'pm'}(Iy,:)); dx = 1./pm;  %% dx at rho points

%% make everything relative to the mount
hm=-h(2,:); [d,Imnt] = min(h(1,:)); x02 = xg(Imnt); %h(1,Imnt)

%% exp.coef
Tcoef = nc.Tcoef(:)*1e-3;  
rho0 = nc.rho0(:);
T0=0;
grav=9.81;

%%-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
%% compute density and stratification
%% this results in N2 having imag values due to negative NN

temp_rho = squeeze(nc{'temp'}(Itime,:,Iy,:));
rho = rho0*(1-Tcoef*(temp_rho-T0)); NN=NaN*ones(size(rho));
NN(2:end-1,:) = -grav/rho0*(rho(3:end,:)-rho(1:end-2,:))./(zz2(3:end,:)-zz2(1:end-2,:));
for ss=2:size(rho,1)-1
    Isel2=find(NN(ss,:)<0); NN(ss,Isel2) = NaN;
end
N2  = sqrt(NN);

%% extrapolate near boundaries
for i=1:size(zz2,2)
    N2([1 end],i) = interp1(zz2(2:end-1,i),N2(2:end-1,i),zz2([1 end],i),'cubic','extrap');
end

dhdx = gradient(hm,xg); % use for criticallity parameter

% figure
% pcolor((xx-x02)/1e3,zz2,real(N2)); shading interp; colorbar
% axis([-50 50 -3000 0])

%% find largest dhdx/alpha
%Isel = find((xg-x02)/1e3>-40 & (xg-x02)/1e3<40);
k=0; alps=[]; rat=[]; N3=[]; II=[];
for i=1:1:length(Isel)
    k=k+1;
    N3(k) = N2(2,Isel(i));
    alps(k) = atan(sqrt(omega^2/(N3(k)^2-omega^2)));
    rat(k)  = dhdx(Isel(i))/tan(alps(k));
    II(k) = i;    
end

rat
% figure
% subplot(2,1,1)
% plot((xg(Isel)-x02)/1e3,hm((Isel)),'k.-')
% axis([-10 10 -1e3 0])
% subplot(2,1,2)
% plot((xg(Isel(II))-x02)/1e3,rat,'k.-')
% %axis([-10 10 -4 4])

%% now determine beams
DX=min(diff(xg)); % horz distance

% counts position of alp and rat
%[f,i1]=max(rat); [f,i2]=min(rat);
%i=1; 

% figure
% plot((xg-x02)/1e3,hm,'k-')
% hold;
% axis([-50 50 -3e3 0])


for i=1:length(rat)
%for i=i1:i2
%for i=1

    %% get quadrants values
    %alp2 = atan2(alps(1),1)
%     AS=1; alp2(AS)=atan2(alps(i),1); SN(AS) = 1;   %1 qd
%     AS=2; alp2(AS)=atan2(alps(i),1); SN(AS) = -1;  %2 qd
%     AS=3; alp2(AS)=atan2(alps(i),1); SN(AS) = -1;  %3 qd
%     AS=4; alp2(AS)=atan2(alps(i),1); SN(AS) = 1;   %4 qd 
    
    AS=1; alp2(AS)=alps(i); SN(AS) = 1;   %1 qd
    AS=2; alp2(AS)=alps(i); SN(AS) = -1;  %2 qd
    AS=3; alp2(AS)=alps(i); SN(AS) = -1;  %3 qd
    AS=4; alp2(AS)=alps(i); SN(AS) = 1;   %4 qd 

    xpos=[]; ypos=[]; xi=[]; yi=[];
    Istart(1:4)=Isel(II(i)); yi(1:4) = hm(Istart); xi(1:4) = xg(Istart);
    xpos(1:4,1) = xi(1:4); ypos(1:4,1) = yi(1:4);
    for j=1:100
        AS=1; alp2(AS) = alp2(AS);         
        AS=2; alp2(AS) = pi - alp2(AS);     
        AS=3; alp2(AS) = alp2(AS) - pi;         
        AS=4; alp2(AS) = - alp2(AS);     

        DZ = tan(alp2)*DX.*SN;
        Istart = Istart+SN;     

        %% new pos
        for AS=1:4
            yi(AS) = yi(AS)+DZ(AS); xi(AS) = xi(AS)+DX*SN(AS); 
            if yi(AS)<hm(Istart(AS)) | yi(AS)>0 ; yi(AS)=NaN; end        
            xpos(AS,j+1) = xi(AS); ypos(AS,j+1) = yi(AS); %log positions

            %% interp new N
            Ni(AS) = interp1(zz2(:,Istart(AS)),N2(:,Istart(AS)),yi(AS),'linear','extrap');
    %        Ni(AS) = interp1(zz2(:,Istart(AS)),N2(:,Istart(AS)),yi(AS));
            alp2(AS) = atan(sqrt(omega^2/(Ni(AS)^2-omega^2)));
%            alp2(AS) = atan2(sqrt(omega^2/(Ni(AS)^2-omega^2)),1);
        end
    end

%     AS=1; plot((xpos(AS,:)-x02)/1e3,ypos(AS,:),'r.-'); 
%     AS=2; plot((xpos(AS,:)-x02)/1e3,ypos(AS,:),'b.-')
%     AS=3; plot((xpos(AS,:)-x02)/1e3,ypos(AS,:),'k.-')
%     AS=4; plot((xpos(AS,:)-x02)/1e3,ypos(AS,:),'g.-')
end
