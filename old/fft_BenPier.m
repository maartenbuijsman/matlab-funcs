%% [X,Xh,fk] = fft_BenPier(t,x)
%% Maarten Buijsman, NIOZ, 25-05-06
%% function fft according Bendat and Piersol (1986)
%% input: vectors t (time) and x (data)
%% note that t=n*dt, n=0,1,2,...N-1, N is length data
%% output: complex fourier components X without filtering and Xh with a Hanning window
%% and frequencies fk=k/(N*dt) n=0,1,2,...N/2, fn = Nyquist frequency (this is the cutoff frequency)
%% first values in X and Xh are the sums of the (smoothed) data

function [X,Xh,fk] = fft_BenPier(t,x);

%% test
% x = [0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0];
% t = [0:length(x)-1];
% t = [0:1:1000]; T=20;
% x = 1*cos(2*pi/T*t) + 2*sin(pi/T*t) ;


N = length(x);
dt = t(2)-t(1);
T = N*dt;
fk = [0:1:floor(N/2)]/T;

%% get fft, same method as matlab fft
%% alternative and faster trick without the loop
%% frequency horizontal and data vertical
k = [0:1:floor(N/2)]; %% over freqs
n = [0:1:N-1];        %% over data points

xm  = ones(length(k),1)*x;  %% data
nm  = ones(length(k),1)*n;  %% data
km  = k'*ones(1,length(x)); %% freqs 

%% raw fft
X  = dt*sum(xm.*exp(-j*2*pi.*km.*nm/N),2); X = X.';
%% with hanning (tapering function)
Xh = dt*sqrt(8/3)*sum(xm.*(1-(cos(pi*nm/N)).^2).*exp(-j*2*pi.*km.*nm/N),2); Xh = Xh.';

%% old loop thing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives the same answ as above
% X = zeros(size(fk)); Xh = X;
% for k=0:1:floor(N/2)      %% over freqs
%     n = [0:1:N-1]; %% over data points
%     X(k+1) = dt*sum(x.*exp(-j*2*pi*k.*n/N)); 
%     
%     %% with hanning (tapering function)
%     Xh(k+1) = dt*sqrt(8/3)*sum(x.*(1-(cos(pi*n/N)).^2).*exp(-j*2*pi*k.*n/N)); 
% end


return


A = abs(xm.*exp(-j*2*pi.*km.*nm/N))
B = abs(exp(-j*2*pi.*km.*nm/N))
figure; surf(A)
figure; surf(B)

A = real(xm.*exp(-j*2*pi.*km.*nm/N))
B = imag(xm.*exp(-j*2*pi.*km.*nm/N))
figure; surf(B)

figure
for i=1:length(fk)
    plot(A(i,:) + B(i,:))
    pause
end

% compare the test with fft from Matlab => perfect 1 on 1!
figure
plot(t,x)


XX = X.*conj(X);

figure
plot(fk,XX)

X2 = fft(x);
XX2 = X2.*conj(X2)

hold
plot(fk,XX2(1:length(fk)),'r.--')


