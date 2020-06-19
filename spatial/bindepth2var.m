function [modelc,obsc,depc] = bindepth2var(model,obs,AA,DD,depb,IO);
%% [modelc,obsc,depc] = bindepth2var(model,obs,AA,depb,
%  MCB, USM, 2020-3-25
%  computes area averaged depth-bin values
%  AA is area; DD is depth
%  IO=1 => display values

% % test
% model = vary(:);
% obs   = varx(:);
% AA    = AA(:);
% % test

%% remove NaN and Inf values
Isel = find(isnan(obs(:)) | isinf(abs(obs(:))) | isnan(model(:)) | isinf(abs(model(:))));
model(Isel) = [];
obs(Isel)   = [];
AA(Isel)    = [];
DD(Isel)    = [];

%% center depths
depc = depb(1:end-1)/2 + depb(2:end)/2;

%% loop ober depth bins
obsc   = [];
modelc = [];
for i=1:length(depb)-1
    Isel = find(DD>depb(i) & DD<=depb(i+1));
    
    obsc(i)   = sum(AA(Isel).*obs(Isel))/sum(AA(Isel));
    modelc(i) = sum(AA(Isel).*model(Isel))/sum(AA(Isel));
end

%% plot if necessary
if IO==1
    figure
    [xx,yy]=stairsc(depb/1e3,depc/1e3,log10(obsc));   plot(xx,yy,'b-','linewidth',1)
    %[xx,yy]=stairsc(depb/1e3,depc/1e3,obsc);   plot(xx,yy,'b-','linewidth',1)
    hold
    [xx,yy]=stairsc(depb/1e3,depc/1e3,log10(modelc)); plot(xx,yy,'r-','linewidth',1)
    %[xx,yy]=stairsc(depb/1e3,depc/1e3,modelc); plot(xx,yy,'r-','linewidth',1)
    legend('alt','model')
end
    
