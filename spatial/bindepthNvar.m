function [DATC,depc] = bindepthNvar(DAT,AA,DD,depb,IO);
%% [DATC,depc] = bindepthNvar(DAT,AA,DD,depb,IO);
%  MCB, USM, 2020-3-25
%  computes area averaged depth-bin values for all variables stored 
%  in 3D array DAT; output in 3D DATC
%  each layer (3rd dim) is a new variable; 2nd dim can be empty
%  AA is area; DD is depth
%  IO=1 => display values

% DAT = varx0; %test

%% remove NaN and Inf values
[ny,nx,nz]=size(DAT);

% locate bad data
ST = [];
Isel = [];
for i=1:nz
    ST(i).vv = DAT(:,:,i);
    Isel = [Isel; find(isnan(ST(i).vv(:)) | isinf(abs(ST(i).vv(:))))];
end

% remove bad data for all variables
for i=1:nz
    ST(i).vv(Isel) = [];
end

AA(Isel)    = [];
DD(Isel)    = [];

%% center depths
depc = depb(1:end-1)/2 + depb(2:end)/2;

%% loop ober depth bins
obsc   = [];
modelc = [];
for i=1:length(depb)-1
    Isel = find(DD>depb(i) & DD<=depb(i+1));
    
    for ii=1:nz
        DATC(ii,i) = sum(AA(Isel).*ST(ii).vv(Isel))/sum(AA(Isel));        
    end
    
end

%% plot if necessary
if IO==1
    figure
    for ii=1:nz
        [xx,yy]=stairsc(depb/1e3,depc/1e3,log10(DATC(ii,:))); plot(xx,yy,'-','linewidth',1)
        holder
    end
end
    
