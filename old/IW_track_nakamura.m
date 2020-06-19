% 
clear all

U0=1; 
T=1; 
t=[0:0.05:2*T];
om=2*pi/T;
U = U0*sin(om*t)+1*U0;
k=2*pi

%figure; plot(t,U)

k*U0/om


figure
for i=1:length(t)
 %i=6
    %% phase speeds of wave
    cplus = -U(i) + om/k;
    cmin  = -U(i) - om/k;
    
    xplus=ones(size(t))*NaN;
    xmin=ones(size(t))*NaN;    
    for j=i+1:length(t)
    
        %% distance of wave advected by current
        xplus(j) = (t(j)-t(i))*cplus + trapz(t(i:j),U(i:j));
        xmin(j) = (t(j)-t(i))*cmin + trapz(t(i:j),U(i:j));        
    end
    subplot(1,2,1)
    plot(xmin,t,'k-')
    if i==1; hold; end

    subplot(1,2,2)
    plot(xplus,t,'k-')
    if i==1; hold; end    
    
end
title(num2str(k*U0/om))