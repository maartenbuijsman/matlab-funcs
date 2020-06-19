function [XC,YC,CNT,IBND] = bincounts2var(model,obs,nbin,IO);
%% [XC,YC,CNT,IBND] = bincounts2var(model,obs,nbin,IO);
%  MCB, USM, 2020-3-26
%  counts the number of variable pairs in one bin
%  IO=1 => display values in figure

%% test
% model = log10(vary(:));
% obs   = log10(varx(:));
% nbin  = 50;
%% test

%% remove NaN and Inf values
Isel = find(isnan(obs(:)) | isinf(abs(obs(:))) | isnan(model(:)) | isinf(abs(model(:))));
model(Isel) = [];
obs(Isel)   = [];

%% bin the data 
% cell faces
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

Imin = max([min(find(sum(CNT,1)>0)) min(find(sum(CNT,2)>0))]);
Imax = min([max(find(sum(CNT,1)>0)) max(find(sum(CNT,2)>0))]);
IBND = [Imin Imax];

% replace zeros with NaNs
CNT(CNT==0) = NaN;


%% Create a figure
if IO==1;
    rcorr = corrcoef_NaN(model(:),obs(:));
    
    figure
    mypcolor(XC,YC,CNT); shading flat
    axis equal
    colorbar
    title(['r = ' num2str(rcorr)])
end


