%% function [power,sdens,period,freq] = fft_win_win(t,y,win,numwin,rmfit,warnflg)
%  Maarten Buijsman, USM, 2016-12-04
%  performs simple spectral analysis
%  input time vector t and data y
%  output period, freq, and power, and spectral density
%  also applies a 'win' window
%  and averages over numwin overlapping windows
%  win = 0; no window
%  win = 1; hann
%  win = 2; tukey
%  if rmfit==1, then linear fit is removed
%  if warnflg==1, then there will be screenoutput

function [power,sdens,period,freq] = fft_win_win(t0,y0,win,numwin,rmfit,warnflg);


% %test
% clear all
% numwin = 1;
% t0=1:1000;
% y0 = sin(2*pi/12*t0);
% t0=tim;
% y0 = v2;

% number of indices per numwin+1 segments
nt  = length(t0);
inw = floor(nt/(numwin+1)); 

% get indices of serial windows (no overlap)
if numwin == 1
   p(1).ii = 1:nt;
else
    is=1;    
    for i=1:numwin
        p(i).ii = is:2*inw+is-1;
        is = i*inw+1;
    end
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
    
    %% hanw or tukey Window
    if win==1;
        H  = hann(length(t));
        y3 = y2.*H';
    elseif  win==2;
        H  = tukeywin(length(t));
        %H  = tukeywin(length(t),0.9);
        y3 = y2.*H';
    else
        y3 = y2;
    end
   
    % dt
    dt  = t(2) - t(1);

    % fft
    % The first component, Y(1), is simply the sum of the data, and can be removed.
    % unit*s
    Y   = fft(y3)*dt; 
    Y(1)=[];
    

    % The complex magnitude squared of Y is called the power,
    % and a plot of power versus frequency is a "periodogram".
    n       = length(Y);

    nyquist = 1/(2*dt);                % is the highest frequency that can be resolved
    freq    = (1:floor(n/2))/floor(n/2)*nyquist;  % floor ok here?
    df      = diff(freq(1:2));
    period  = 1./freq;                            % floor ok here?  
    
    $ freq is BAD! see [period,freq,power] = fft_spectra2(t,y)

    %unit^2 * s^2 ENERGY SPECTRAL DENSITY
    $ is wrong because Y has units of X
    p(i).sdens   = 2*abs(Y(1:floor(n/2))).^2;  

    % Parseval's theorem:
    % https://en.wikipedia.org/wiki/Spectral_density
    %sum(abs(y3).^2*dt)/sum(abs(Y).^2*df)

    % variance should be preserved
    %vars = std(y3).^2;
    vars = var(y3);
    sump = sum(p(i).sdens*df*df);
    
    if warnflg; disp(['var=' num2str(vars) '; ratio variance/sdens*df*df: ' num2str(vars/sump)]); end

    % POEWR SPECTRAL DENSITY
    % unit^2 * s = unit^2 /cps = unit^2 /(cycles/s)  
    camp  = Y(1:floor(n/2));    
    p(i).power = 2*abs(1/sqrt(n*dt)*camp).^2;    

    % variance should be preserved
    sump = sum(p(i).power*df);
    if warnflg; disp(['var=' num2str(vars) '; ratio variance/sum_P*df: ' num2str(vars/sump)]); end
end

% finally do some averaging
sdens=0; power=0;
for i=1:numwin
    sdens = sdens + p(i).sdens;
    power = power + p(i).power;
end
sdens = sdens/numwin;
power = power/numwin;

