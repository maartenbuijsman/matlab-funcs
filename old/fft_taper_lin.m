%% function [period,freq,power,camp] = fft_taper_lin(t,y,apct,rmfit)
%  Maarten Buijsman, USM, 2016-03-15
%  performs simple spectral analysis
%  input:  equidistant time vector t and data y
%          apct is percentage of total time series that needs cosine tapering; 30 is a good value
%          5 gives weird lobe oscillations
%          if apct=0, no tapering occurs 
%          when rmfit==1, the linear fit is removed from the time series
%  output: period, freq [cp unit time],
%          complex amplitudes camp
%          and power [unit_y^2 unit_time]
%
%  based on http://www.mathworks.com/products/demos/matlab/sunspots/sunspots.html
%  see also help_understand_fft.m, where this function is check
%  note that the longer the time series, the better the resolution

function [period,freq,power,camp] = fft_taper_lin(t,y1,apct,rmfit);

% % test
% apct=5;
% rmfit = 1;
% t=0:1:700;
% y1 = cos(2*pi/12.5*t) + 0.25*cos(2*pi/4*t);

nx = length(t);  % length of time series
ii = 1:nx;       % indices

%% 1) remove linear trend
if rmfit==1
    cf    = polyfit(ii,y1,1);
    y_fit = polyval(cf,ii);
    y2 = y1 - y_fit;    
else
    y2 = y1;
end

%figure;  plot(ii,y1,'k-',ii,y2,'r--') %compare original and modified time series

%% 2) cosine taper
% 1/2*(1-cos(pi/a*t)) for t<a
% 1/2*(1-cos(pi/a*(t(end)-t))) for t > t-a
y3 = y2;

a      = apct/100*nx;  % in time series index
I1     = find(ii<a);  % taper start
y3(I1) = 1/2*(1-cos(pi/a*ii(I1))).*y2(I1);

I2     = find(ii>ii(end)-a);   % taper end
y3(I2) = 1/2*(1-cos(pi/a*(ii(end)-ii(I2)))).*y2(I2);

%figure;  plot(1/2*(1-cos(pi/a*(ii(end)-ii(I2))))) %check
%figure;  plot(ii,y1,'k-',ii,y3,'r--') %compare original and modified time series

%% now extract the compex amplitudes 
%  which include the amplitudes and phases
Y = fft(y3);

% The first component of Y, Y(1), is simply the sum of the data, and can be removed.
Y(1)=[];

% The complex magnitude squared of Y is called the power,
% and a plot of power versus frequency is a "periodogram".

n       = length(Y);
power   = abs(Y(1:floor(n/2))).^2;
camp    = Y(1:floor(n/2));

dt      = t(2)-t(1);
nyquist = 1/(2*dt);                % is the highest frequency that can be resolved
freq    = (1:floor(n/2))/floor(n/2)*nyquist;  % floor ok here?
period  = 1./freq;                            % floor ok here?  

return

%test

%plot
figure; 
loglog(freq*24,power,'k-')
set(gca,'xtick',[0.5 1 2 4 6])

[X,Xh,fk] = fft_BenPier(t,y1);

hold
loglog(fk*24,Xh.*conj(Xh),'r-');
loglog(fk*24,X.*conj(X),'b-');

% hann window
H  = hann(length(y1));  % samples of a Hann window 

%figure; plot(t,y1.*H')
Y2 = fft(y1.*H');
Y2 = Y2(1:length(freq)+1)

loglog(fk*24,Y2.*conj(Y2),'g-');

