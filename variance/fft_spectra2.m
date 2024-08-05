%% function [period,freq,power] = fft_spectra2(t,y,tukey,numwin,linfit)
%  Maarten Buijsman, USM, 2024-7-23
%  performs simple spectral analysis
%  input
%     time vector t1 and data y1
%     tukey and numwin 
%  output 
%     period, freq, and power [y_unit^2/(cycle/t_unit)]
%     e.g. y_unit = m and t_unit = day --> [m^2/cpd]
%
%  if tukey = 1, a hanning window is applied 
%  if tukey = 0, a boxcar window is applied (all ones; no tukey window)
%  if tukey < 1, a tukey window is applied with R = tukey
%
%  numwin sets the number of 50% overlapping windows 
%  if numwin = 1, the entire even length of the time series is used
%  the time series mean is removed
%  if linfit = 1, the linear fit is removed
%
%  In the code Parseval's theoreum is satisfied: 
%     sum(abs(xi).^2)*dt = sum(abs(Xi).^2)*df 
%     see https://www.mathworks.com/matlabcentral/answers/15770-scaling-the-fft-and-the-ifft 

function [period,freq,power] = fft_spectra2(t1,y1,tukey,numwin,linfit);

%% test
% clear all
% numwin = 3
% tukey=0;
% T1 = 12; T2=24; T3=48;
% t1 = 1:1447;
% y1 = 1*cos(2*pi/T1*t1) + 0.5*cos(2*pi/T2*t1) + 0.25*cos(2*pi/T3*t1);
% figure; plot(t1,y1)

% inw is half length of windowed time series 
% thus the length of the time series is always even
nt1 = length(t1);
inw = floor(nt1/(numwin+1)); 

% store indices of each window in p(i)
% these are 50% overlapping windows 
is=1;
for i=1:numwin
    p(i).ii = is:2*inw+is-1;
    is = i*inw+1;
end

% do fft for each window
for i=1:numwin
    
    y = y1(p(i).ii);
    t = t1(p(i).ii);  
    nt = length(t);

    % remove linear fit (this is a slow routine ...)
    if linfit
        cf    = polyfit(t-mean(t),y,1);
        y_fit = polyval(cf,t-mean(t));
        y = y - y_fit;  
    end

    % fast!
    y = y-mean(y);

    % use a Tukey window
    H = tukeywin(length(t),tukey);
    y = y.*H';

    % freq: 1st value is 1/(nt*dt), last value is fn = 1/(2*dt)
    dt    = t(2) - t(1);
    df    = 1/(dt*nt);
    freq  = 1/dt*(1:(nt/2))/nt;

    % Y(1) is the sum and is omitted
    % Y(2) = conjugate(Y(end))
    % index nt/2+1 is shared between left and right side of spectrum
    % real numbers for these sides are the same 
    % imaginary numbers for these sides are the complex conjugates
    Y  = fft(y)*dt;           % y_unit*s

    % remove the sum value, now nt/2 is the shared value
    Y(1) = [];                

    % Parseval's theorem; ratio between integrated energy and integrated spectral denisty  (should be 1)
    %sum(y.^2*dt)/sum(abs(Y).^2*df)

    P2 = abs(Y).^2;           % energy y_unit^2*s^2
    P1 = 2*P2(1:nt/2);        % select left side and double the power (= folding the spectrum)
    powers(i,:) = P1*df;      % y_unit^2*s^2 * 1/s = y_unit^2*s = y_unit^2*1/Hz      (Hz = 1/s)

    % Parseval's theorem ratio between integrated energy and power (should be 1)
    %sum(y.^2*dt)/sum(powers)

    % omit first value
    period  = 1./freq;  
end

% average the power over the number of windows
power = mean(powers,1);

% test
%figure; semilogy(freq*24,power)



