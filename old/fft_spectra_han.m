%% function [period,freq,power] = fft_spectra_han(t,y)
%% Maarten Buijsman, NIOZ, 2010-05-07
%% performs simple spectral analysis
%% input time vector t and data y
%% output period, freq, and power [unit_y^2 unit_time]
%% based on http://www.mathworks.com/products/demos/matlab/sunspots/sunspots.html
%% see also help_understand_fft.m, where this function is check
%% note that the longer the time series, the better the resolution
%% also applies a hanning window

function [period,freq,power,powerH] = fft_spectra_han(t,y);

Y = fft(y);

%% The first component of Y, Y(1), is simply the sum of the data, and can be removed.
Y(1)=[];

%% The complex magnitude squared of Y is called the power, 
%% and a plot of power versus frequency is a "periodogram".

n=length(Y);
power = abs(Y(1:floor(n/2))).^2;
dt = mean(diff(t));
nyquist = 1/(2*dt);  %% highest frequency that can be resolved: 1/(2(dt)); dt=1 yr
freq = (1:n/2)/(n/2)*nyquist;
period=1./freq;

%figure; plot(freq,power); xlabel('cycles/day'); title('Periodogram');
%figure; plot(period,power); ylabel('Power'); xlabel('Period (days/Cycle)');


%% does the hanning trick
H  = hann(length(y));  % samples of a Hann window 
Y2 = fft(y.*H');

%% The first component of Y, Y(1), is simply the sum of the data, and can be removed.
Y2(1)=[];
powerH = abs(Y2(1:floor(n/2))).^2;
