function [k,L,C,Cg,Ce,Weig,Ueig] = sturm_liouville_noneqDZ_norm(zf,N2,f,om,nonhyd)
%% [k,L,C,Cg,Ce,Weig,Ueig] = sturm_liouville_noneqDZ_norm(zf,N2,f,om,nonhyd) 
%
% Maarten Buijsman & Oladeji Siyanbola, USM, 2024-12-28
% computes eigenvalues using Ashok & Bhaduria (2009) numerical scheme
%
% input:
% N2 [rad^2/s^2] Brunt Vaisala frequency-squared at layer faces, surface to bottom
% zf [m] layer faces' depth vector (approaches -inf for increasing depth)
% f [rad/s] coriolios frequence; om [rad/s] internal wave frequency
% nonhyd =-1 then solve nondyrostatic SL eq. 
%
% output:
% k [rad/m], L [m],C [m/s],Cg [m/s], Ce [m/s] wavenumber, wavelength, phase speed, group speed, eigen speed
% Weig [m], Ueig [-] vert. (@ layer interfaces) and hor. (@ cell centers) eigen functions 

dz = -diff(zf); % layer thickness vector (top to bottom)

% Compute parameters
H = nansum(dz);
N = numel(dz);

if nonhyd
    %% non-hydrostatic part: N^2-omega^2     
    NN = N2-om^2; 
    Isel=find(NN<0); 
    NN(NN<0) = 1e-12;
else
    %% hydrostatic part
    NN = N2;     
end
B = diag(-NN(2:end-1));

% Main solver
A = [];
A(1,1) = -2*1/(dz(1)*dz(2));
A(1,2) = 2*1/(dz(2)*(dz(1) + dz(2)));
for i = 2:N-2
    A(i,-1+i) =  2*1/(dz(i)*(dz(i) + dz(i+1)));
    A(i, 0+i) = -2*1/(dz(i)*dz(i+1));
    A(i, 1+i) =  2*1/(dz(i+1)*(dz(i) + dz(i+1)));
end
i = i+1;
A(i,-1+i) =  2*1/(dz(i)*(dz(i) + dz(i+1)));
A(i, 0+i) = -2*1/(dz(i)*dz(i+1));

ll = size(A,1);

% Solve the EVP
[W1,invCe2] = eig(A,B);
Ce2 = 1./(diag(invCe2));
Ce = sqrt(Ce2); [Ce,Is]=sort(Ce,'descend');
k = abs(sqrt((om^2-f^2))./Ce);
W1 = W1(:,Is);
C  = om./k;
L  = 2*pi./k;
Cg = Ce.^2.*k/om;
%% Get W structure functions
W2 = [zeros(1,ll); W1; zeros(1,ll)];    % W at cell faces

%% Get U structure
dW2 = W2(2:end,:)-W2(1:end-1,:);
dzu = repmat(dz,[1,size(dz,1)-1]);
Ueig1 = dW2./dzu;
AA = repmat(sum(Ueig1.^2.*dzu,1)./H,[N 1]).^(1/2);
AA(AA==0) = Inf;
Ueig2 = Ueig1./AA;
Ueig2(:,Ueig2(N,:)<0) = -Ueig2(:,Ueig2(N,:)<0);

% Prepare for Weig & Ueig
Ueig = Ueig2;
Weig = W2;
