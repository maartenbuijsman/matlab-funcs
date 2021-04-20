%% function [period,freq,power] = fft_spectra2(t,y)
%  Maarten Buijsman, USM, 2020-8-25
%  performs simple spectral analysis
%  input time vector t and data y
%  output period, freq, and power [y_unit^2*1/Hz]
%  added window and trend removal
%  fixed frequencies <= works for even length time series
%  Made sure Parseval's theoreum is used: sum(abs(xi).^2)*dt = sum(abs(Xi).^2)*df 
%  see https://www.mathworks.com/matlabcentral/answers/15770-scaling-the-fft-and-the-ifft 

function [period,freq,power] = fft_spectra2(t,y);

% % %% test
% clear all
% T1 = 12; T2=24; T3=35;
% t = 1:720;
% y = 1*cos(2*pi/T1*t) + 0.5*cos(2*pi/T2*t) + 0.25*cos(2*pi/T3*t);
% figure; plot(t,y)

% make sure the time series is even
n1 = length(t);
if rem(n1,2) ~= 0
    t = t(1:end-1);
    y = y(1:end-1);
    n1 = length(t);
end

%whos t y

% % remove linear fit (this is a slow routine ...)
% cf    = polyfit(t,y,1);
% y_fit = polyval(cf,t);
% y = y - y_fit;  

% fast!
%y = y-mean(y);

% % use a Tukey window
% H = tukeywin(length(t));
% y = y.*H';

% freq: 1st value is 1/(n1*dt), last value is fn = 1/(2*dt)
dt    = t(2) - t(1);
df    = 1/(dt*n1);
freq  = 1/dt*(0:(n1/2))/n1;

% Y(1) relates to the mean
% Y(2) = conjugate(Y(end))
% index n1/2+1 is shared between left and right side of spectrum
% real numbers for these sides are the same 
% imaginary numbers for these sides are the complex conjugates
Y  = fft(y)*dt;               % y_unit*s

% Parseval's theorem; ratio between integrated energy and integrated spectral denisty  (should be 1)
% sum(y.^2*dt)/sum(abs(Y).^2*df)

P2 = abs(Y).^2;               % energy y_unit^2*s^2 = y_unit^2*1/Hz^2  
P1 = P2(1:n1/2+1);            % select left side (1:n+1), keep first value, which relates to the mean 
P1(2:end-1) = 2*P1(2:end-1);  % double the power (= folding the spectrum)
power = P1*df;                % y_unit^2*s^2 * 1/s = y_unit^2*s = y_unit^2*1/Hz      (Hz = 1/s)

% check; these values should be the same!!
% sum(power)
% sum(abs(y).^2*dt)

% omit first value
power = power(2:end);
freq  = freq(2:end);
period  = 1./freq;                        

% test
%figure
%semilogy(period,power)



