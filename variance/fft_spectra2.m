%% function [period,freq,power] = fft_spectra2(t,y,tukey)
%  Maarten Buijsman, USM, 2021-4-20
%  performs simple spectral analysis
%  input time vector t and data y
%  output period, freq, and power [y_unit^2*1/Hz]
%  if tukey = 1, a tukeywindow is applied 
%  mean is removed
%  fixed frequencies <= works for even length time series
%  Made sure Parseval's theoreum is used: sum(abs(xi).^2)*dt = sum(abs(Xi).^2)*df 
%  see https://www.mathworks.com/matlabcentral/answers/15770-scaling-the-fft-and-the-ifft 

function [period,freq,power] = fft_spectra2(t,y,tukey);

% %% test
% clear all
% tukey=0;
% T1 = 12; T2=24; T3=48;
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

% % remove linear fit (this is a slow routine ...)
% cf    = polyfit(t,y,1);
% y_fit = polyval(cf,t);
% y = y - y_fit;  

% fast!
y = y-mean(y);

% use a Tukey window
if tukey
    H = tukeywin(length(t));
    y = y.*H';
end

% freq: 1st value is 1/(n1*dt), last value is fn = 1/(2*dt)
dt    = t(2) - t(1);
df    = 1/(dt*n1);
freq  = 1/dt*(1:(n1/2))/n1;

% Y(1) is the sum and is omitted
% Y(2) = conjugate(Y(end))
% index n1/2+1 is shared between left and right side of spectrum
% real numbers for these sides are the same 
% imaginary numbers for these sides are the complex conjugates
Y  = fft(y)*dt;           % y_unit*s

% remove the sum value, now n1/2 is the shared value
Y(1) = [];                

% Parseval's theorem; ratio between integrated energy and integrated spectral denisty  (should be 1)
%sum(y.^2*dt)/sum(abs(Y).^2*df)

P2 = abs(Y).^2;           % energy y_unit^2*s^2
P1 = 2*P2(1:n1/2);        % select left side and double the power (= folding the spectrum)
power = P1*df;            % y_unit^2*s^2 * 1/s = y_unit^2*s = y_unit^2*1/Hz      (Hz = 1/s)

% Parseval's theorem ratio between integrated energy and power (should be 1)
%sum(y.^2*dt)/sum(power)

% omit first value
period  = 1./freq;                        

% test
%figure
%semilogy(freq*24,power)



