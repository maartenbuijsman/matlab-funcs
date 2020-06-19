%% function [Ueig,Weig,Beig,Peig]=bfq_eigenvectors(rhoNil,fcor,omfq,Ieig,zw,zr,Nbfq,Cf,Lf,C,L,gam,V,plotfig);         
%% compute eigenvecors
%% mcb, gfdl, 2011-05-05
%% firs use  [Cf,Lf,C,L,gam,V] = sturm_liouville_nonhyd(omfq2,zw,Nbfq,fcor); 

function  [Ueig,Weig,Beig,Peig]=bfq_eigenvectors(rhoNil,fcor,omfq,Ieig,zw,zr,Nbfq,Cf,Lf,C,L,gam,V,plotfig);        

[C2,Ia] = sort(C,'descend');                    % sort C1 in descending order
Cf = Cf(Ia);
V2 = V(:,Ia); L = L(Ia); Lf = Lf(Ia);    % re-arrange
if     fcor~=0;   kn = 2*pi./Lf(1:Ieig);    % used for Normalizing WITH f
elseif fcor==0;   kn = 2*pi./L(1:Ieig);     % used for Normalizing NO f
end

Weig=[]; Weig(2:size(V2,1)+1,:) = V2(:,1:Ieig); % select Ieig number of eigenmodes
Weig(end+1,:) = 0;

%% W eigenvectors -----------------------------------------------------------------------
%% test to check ORTHOGONALITY
%% indeed for I1 is not I2 => eigsum~0
%% for  I1=I2, the function needs to be normalized!!! to 1
%if NONhyd == false
%            R     = (NS(:,Ix(JJ)).^2)/(omegasel^2-f^2);  % Wzz + k2"R"W =   0
%elseif NONhyd == true
%    R     = (NS(:,Ix(JJ)).^2-omegasel^2)/(omegasel^2-f^2);  % Wzz + k2"R"W =   0
%end
%R     = (Nbfq.^2)/(omfq^2-fcor^2);  % hydrostatic
R     = (Nbfq.^2-omfq^2)/(omfq^2-fcor^2);  % NON-hydrostatic
knorm = sqrt(1/( trapz(zw,R.*Weig(:,1).*Weig(:,1))/abs(zw(1)) ) );     % normalizing factor, depth-averaged
I1 = 6; I2 = 7;                                        % test
eigsum = cumtrapz(zw,R.*Weig(:,I1).*Weig(:,I2)*knorm^2);
eigtot = trapz(zw,R.*Weig(:,I1).*Weig(:,I2)*knorm^2);
Weig   = Weig*knorm;                                   % normalize Weig
Weig2  = Weig;                                         % at w-points; for Ueig and Beig
Weig   = (Weig(2:end,:)+Weig(1:end-1,:))/2;            % get W at rho points

%% keep on same side
Isel = find(Weig(1,:)<0);
Weig(:,Isel) = -Weig(:,Isel);
if plotfig;  figure; plot(Weig(:,[1:3]),zr); end

%% U eigenvectors -----------------------------------------------------------------------
%% if U = dW/dz *1/kn is used then no norm factor needs to be computed
%% else
DZ   = diff((zw*ones(1,size(Weig,2))),[],1);
Ueig = diff(Weig2,[],1)./DZ.*1./(ones(size(DZ,1),1)*kn');  % qn is constant

%% normalize the eigen vector by integrating the square
%% and normalizing over the waterdepth (ASIAEX Ying-Jang 2004)
%% do not need R in int R*U*U=0/1 ?? why not ??
thisisone = sum(Ueig.*Ueig.*DZ,1)/abs(zw(1));   % check!
thisiszero= sum(Ueig(:,1).*Ueig(:,2).*DZ(:,1));         % check!

%% keep on same side
Isel = find(Ueig(1,:)<0);
Ueig(:,Isel) = -Ueig(:,Isel);
if plotfig; figure; plot(Ueig(:,[1:3]),zr); end

%% Buoyancy eigenvectors -----------------------------------------------------------------------
Beig   = -1/omfq*repmat(Nbfq,[1 Ieig]).^2.*Weig2;

%         %% check, why normalize by 1/N^2 to obtain 1/zero => to get W, dude!
%         thisisone  = trapz(zw,1./NS(:,Ix(JJ)).^2.*Beig(:,1).*Beig(:,1))/abs(hm(Ix(JJ))); %% nearly zero.... 
%         thisiszero = trapz(zw,1./NS(:,Ix(JJ)).^2.*Beig(:,1).*Beig(:,2))/abs(hm(Ix(JJ))); %% nearly zero....         
%         eigsum     = cumtrapz(zw,1./NS(:,Ix(JJ)).^2.*Beig(:,1).*Beig(:,2))/abs(hm(Ix(JJ)));
%         %figure; plot(eigsum,zw)

Beig   = Beig(2:end,:)/2+Beig(1:end-1,:)/2; %at rho pints

%% keep on same side
Isel = find(Beig(1,:)<0);
Beig(:,Isel) = -Beig(:,Isel);
if plotfig;  figure; plot(Beig(:,[1:3]),zr); end

%% P eigenvectors -------------------------------------------
Peig = rhoNil*(omfq.^2 - fcor^2)./(omfq*(ones(size(DZ,1),1)*kn').^2).*diff(Weig2,[],1)./DZ;            

%% make sure the base is positive => otherwise phase inconsistencies!!     
isel = find(Peig(1,:)<0); 
Peig(:,isel) = -Peig(:,isel);

if plotfig;  figure; plot(Peig(:,[1:3]),zr); end



return
