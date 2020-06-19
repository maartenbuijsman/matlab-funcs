%% function [ZZCa,DZCW,XXC,ZZC] = gridZC(XG,XC,RC,DZC,DRF);
%% MCB, GFDL, 2011-01-07
%% computes new MITgcm grid ZZCa and others for pcolorplot
%% input from ana* or dissip* routines

function [ZZCa,DZCW,XXC,ZZC] = gridZC(XG,XC,RC,DZC,DRF);

XXC = repmat(XC,[length(RC) 1]);
ZZC = repmat(RC,[1 length(XC)]);  

%% to compute gradient of W
DZCW = DZC;
for i=1:size(DZC,1)
    I0 = find(DZC(i,:)==0);
    DZCW(i,I0) = 1;
end

DZa1 = DZCW; DZa2 = DZCW;
for i=1:size(DZCW,1)
    I0 = find(DZCW(i,:)~=1);    DZa1(i,I0) = 0;
    I0 = find(DZCW(i,:)==1);    DZa2(i,I0) = 0;    
    I0 = find(DZCW(i,:)>1);     DZa2(i,I0) = 1;        
end

ZZ1 = [zeros(size(XC)); cumsum(DZC(end:-1:1,:),1)]; %% at faces, adjusted
ss = cumsum(DRF(end:-1:1));
ZZ2 = [zeros(size(XC)); repmat(ss,[1 length(XG)])]; %% at faces, not adjusted

%% reverse
ZZ1b = ZZ1(end:-1:1,:);
ZZ2b = ZZ2(end:-1:1,:);

%% average
ZZ1bav = (ZZ1b(1:end-1,:)+ZZ1b(2:end,:))/2;
ZZ2bav = (ZZ2b(1:end-1,:)+ZZ2b(2:end,:))/2;

%% combine
ZZCa = -ZZ1bav.*DZa2 -ZZ2bav.*DZa1; %% at centers

