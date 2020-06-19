%% function [X1,X2,Xh1,Xh2,fk,y1,y2,t,N,Ns,dt] = crosscor_BenPier(t,y1,y2);
%% Maarten Buijsman, NIOZ, 30-05-06
%% calculates cross correlation R12 and weighted cross correlation rho12 for two variables
%% according to Bendat and Piersol (1986) p.406
%% INPUT: vectors y1 and y2 (data), and percentage of data used proc_N for lagging 
%% OUTPUT: auto correlation R11 and R22 and cross correlation R12 and weighed rho12 (-1<=rho12<=1),
%% and the lag 'time' series r = 0:1:m, m is number of lags (note that tlag = dt*lag)
%% also note that the 'lag' indicates that y2 leads y1

function [R11,R22,R12,rho12,lag] = crosscor_BenPier(y1,y2,proc_N);

% %% test antiphase, at lag = 5, they are in phase
% t=0:99;
% y1 = cos(2*pi/10*t);    %% behind in phase
% y2 = cos(2*pi/10*t-pi); %% ahead in phase (leading)

%proc_N = 70;
N = length(y1); 
m = round(N*proc_N/100); %% number of lags

R12 = ones(1,length(m+1)); lag = R12;
for r=0:m
    n=1:N-r;
    R12(r+1) = 1/(N-r)*sum(y1(n).*y2(n+r));
    lag(r+1) = r;
    %plot(n,y1(n)/std(y1),'r',n,y2(n+r)/std(y2),'b'); pause
end

R11 = 1/N*sum(y1.*y1);
R22 = 1/N*sum(y2.*y2);
rho12 = R12/(sqrt(R11*R22)); %% normalized cross correlation

% figure; plot(lag,rho12)