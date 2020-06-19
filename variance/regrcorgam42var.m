function [Areg,rcorr,gamma] = regrcorgam42var(model,obs,AA,IO);
%% [Areg,rcorr,gamma] = regrcorgam42var(model,obs,AA,IO);
%  MCB, USM, 2020-3-25
%  computes three statistics for comparing model with observations:
%  regression, ratio gamma, correlation
%  (1) linear (2) log10
%  AA is area for the gamma calculation
%  IO=1 => display values
%  previously named as stats42var.m

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


%% regression
% MODEL = OBS*A  PAPER!!!
% A = OBS\MODEL
% A ~ MODEL/OBS
Areg(1) = obs(:)\model(:);
Areg(2) = log10(obs(:))\log10(model(:));


%% ratio gamma model/obs
gamma(1) =  sum(model(:).*AA(:))/sum(AA(:)) ...
      / (sum(obs(:)  .*AA(:))/sum(AA(:))) ;
gamma(2) =  sum(log10(model(:)).*AA(:))/sum(AA(:)) ...
      / (sum(log10(obs(:))  .*AA(:))/sum(AA(:))) ;
  

%% correlation
rcorr(1) = corrcoef_NaN(model(:),obs(:));
rcorr(2) = corrcoef_NaN(log10(model(:)),log10(obs(:)));


%% display results 
if IO
    disp(['reg A is ' num2str(Areg)]);
    disp(['gamma is ' num2str(gamma)]);
    disp(['rcorr is ' num2str(rcorr)]);
end

return

%% why is the below here?
%% bin the data 
nbin = 100;
bnds = linspace(min([model; obs]),max([model; obs]),nbin+1);
[X,Y] = meshgrid(bnds,bnds);

% center coordinates
xc = bnds(1:end-1)/2 + bnds(2:end)/2;
yc = bnds(1:end-1)/2 + bnds(2:end)/2;

[XC,YC] = meshgrid(xc,yc);

CNT = XC*0;
for i=1:nbin
    for j=1:nbin    
        xp = [X(j,i) X(j,i+1) X(j+1,i+1) X(j+1,i) X(j,i)];
        yp = [Y(j,i) Y(j,i+1) Y(j+1,i+1) Y(j+1,i) Y(j,i)];        
        IN = inpolygon(obs,model,xp,yp);
        CNT(j,i) = length(find(IN==1));
    end
end



