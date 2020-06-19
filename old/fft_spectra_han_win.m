%% function [period,freq,power] = fft_spectra_han_win(t,y,numwin)
%  Maarten Buijsman, NRL, 2014-05-29
%  performs simple spectral analysis
%  input time vector t and data y
%  output period, freq, and power [unit_y^2 unit_time]
%  based on http://www.mathworks.com/products/demos/matlab/sunspots/sunspots.html
%  see also help_understand_fft.m, where this function is check
%  note that the longer the time series, the better the resolution
%  also applies a hanning window
%  and averages over numwin overlapping windows

function [power,powerH,period,freq] = fft_spectra_han_win(t0,y0,numwin);

% %test
% clear all
% numwin = 5;
% t1=1:1000;
% y1 = sin(2*pi/12*t1);

% number of indices per numwin+1 segments
nt  = length(t0);
inw = floor(nt/(numwin+1)); 

% get indices
is=1;
for i=1:numwin
    p(i).ii = is:2*inw+is-1;
    is = i*inw+1;
end

for i=1:numwin
    
    y1 = y0(p(i).ii);
    t  = t0(p(i).ii);  

    % remove linear trend
    if rmfit==1
        cf    = polyfit(t,y1,1);
        y_fit = polyval(cf,t);
        y2 = y1 - y_fit;    
    else
        y2 = y1;
    end
    
    %% Hanning Window
    if hanw==1
        H = hann(length(tim));
        y3 = y2.*H';
    else
        y3 = y2;
    end    
   
    % dt
    dt  = t(2) - t(1);

    % fft
    % The first component, Y(1), is simply the sum of the data, and can be removed.
    Y   = fft(y3)*dt; 
    Y(1)=[];

    % The complex magnitude squared of Y is called the power,
    % and a plot of power versus frequency is a "periodogram".
    n       = length(Y);

    nyquist = 1/(2*dt);                % is the highest frequency that can be resolved
    freq    = (1:floor(n/2))/floor(n/2)*nyquist;  % floor ok here?
    df      = diff(freq(1:2));
    period  = 1./freq;                            % floor ok here?  

    %m^2/s^2*s^2 SPECTRAL DENSITY
    p(i).sdens   = 2*abs(Y(1:floor(n/2))).^2;  

    % Parseval's theorem:
    % https://en.wikipedia.org/wiki/Spectral_density
    %sum(abs(y3).^2*dt)/sum(abs(Y).^2*df)

    % variance should be preserved
%     var = std(y2).^2;
%     sump = sum(sdens*df*df);
%     var/sump

    % power spectral density
    % m2/s2 1/cps
    camp  = Y(1:floor(n/2));    
    p(i).power = 2*abs(1/sqrt(n*dt)*camp).^2;    

    % variance should be preserved
%     sump = sum(p(i).power*df);
%     var/sump
end

% finally do some averaging
sdens=0; power=0;
for i=1:numwin
    sdens = sdens + p(i).powerH;
    power = power + p(i).power;
end
sdens = sdens/numwin;
power = power/numwin;

period=p(1).period;
freq=p(1).freq;




