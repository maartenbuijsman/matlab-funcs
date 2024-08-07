%% function [period,freq,power] = fft_spectra2(t,y,tukey,numwin,linfit,prewhit)
%  Maarten Buijsman, USM, 2024-8-7
%  performs simple spectral analysis
%  input
%     time vector t1 and data y1
%     tukey, numwin, linfit, prewhit 
%  output 
%     period, freq, and power [y_unit^2/(cycle/t_unit)]
%     e.g. y_unit = m and t_unit = day --> [m^2/cpd]
%
%  if tukey = 1, a hanning window is applied 
%  if tukey = 0, a boxcar window is applied (all ones; no tukey window)
%  if tukey < 1, a tukey window is applied with R = tukey
%
%  numwin sets the number of 50% overlapping windows 
%  if numwin  = 1, the entire even length of the time series is used
%
%  if linfit  = 1, the linear fit is removed
%
%  if prewhit = 1, fft is done over dy/dt, P is normalized by freq^2 to re-red the spectrum 
%  prewhitening reduces the adulterating effect of low frquencies
%  not sure if freq^2 normalization is the default or somewhat arbitrary (maybe for freq^-2 slopes?)  
%
%  the time series mean is always removed
%
%  In the code Parseval's theoreum is satisfied: 
%     sum(abs(xi).^2)*dt = sum(abs(Xi).^2)*df 
%     see https://www.mathworks.com/matlabcentral/answers/15770-scaling-the-fft-and-the-ifft 

function [period,freq,power] = fft_spectra2(t1,y1,tukey,numwin,linfit,prewhit);

%% test
% clear all
% numwin = 3
% tukey=0;
% T1 = 12; T2=24; T3=48;
% t1 = 1:1447;
% y1 = 1*cos(2*pi/T1*t1) + 0.5*cos(2*pi/T2*t1) + 0.25*cos(2*pi/T3*t1);
% figure; plot(t1,y1)

% make sure time series length is even number!
if rem(length(t1),2)~=0
    t1(end) = [];
    y1(end) = [];    
end

% dt is independent of number of windows
dt    = t1(2) - t1(1);

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
powers = [];
for i=1:numwin
    
    y  = y1(p(i).ii);
    nt = length(y);
    
    % remove linear trend
    if linfit
        y = detrend(y,1);  
    end

    % fast!
    y = y-mean(y);

    % use a Tukey window
    H = tukeywin(nt,tukey);
    y = y.*H';

    % pre-whitening; this may improve high frequency signals
    % http://pordlabs.ucsd.edu/sgille/sio221c/lectures/spectra_windowing.pdf
    % compute the derivative
    if prewhit
        yd = y;          % for plotting purposes
        y  = diff(y)/dt;   
        
        % make sure time series length is even number!
        if rem(length(y),2)~=0
            y(end) = [];    
        end
        nt = length(y);        
 
%         figure
%         plot(y,'r-')
%         hold on
%         plot(yd,'g-')
    end

    % freq: 1st value is 1/(nt*dt), last value is fn = 1/(2*dt)
    df     = 1/(dt*nt);
    freq   = 1/dt*(1:(nt/2))/nt;
    period = 1./freq;  
    
    % Y(1) is the sum and is omitted
    % Y(2) = conjugate(Y(end))
    % index nt/2+1 is shared between left and right side of spectrum
    % real numbers for these sides are the same 
    % imaginary numbers for these sides are the complex conjugates
    Y  = fft(y)*dt;           % y_unit*t_unit

    % remove the sum value, now nt/2 is the shared value
    Y(1) = [];                

    % Parseval's theorem; ratio between integrated energy and integrated spectral denisty  (should be 1)
    %sum(y.^2*dt)/sum(abs(Y).^2*df)

    P2 = abs(Y).^2;               % energy y_unit^2*t_unit^2
    P1 = 2*P2(1:nt/2);            % select left side and double the power (= folding the spectrum)
    
    if prewhit
        % normalize power by frequency^2 to rered the spectrum   
        powers(i,:) = P1*df./freq.^2; 
    else
        % y_unit^2*t_unit^2 * 1/t_unit = y_unit^2*t_unit = y_unit^2*1/Hz (Hz = 1/t_unit)        
        powers(i,:) = P1*df;          
    end

    % Parseval's theorem ratio between integrated energy and power (should be 1)
    % note that re-redding is excluded 
    % sum(y.^2*dt)/sum(P1*df)  % this is OK

    % NOTE that including re-reddening Parseval's theorem does not hold anymore ......
    % sum(yd.^2*dt)/sum(P1*df./freq.^2) 
end

% average the power over the number of windows
power = mean(powers,1);

% % test
% figure; 
% loglog(freq,power,'g-')
% hold on
% loglog([2 2],[mean(power) mean(power)],'r*')
% loglog([4 4],[mean(power) mean(power)],'g*')
% loglog([6 6],[mean(power) mean(power)],'c*')
% loglog([8 8],[mean(power) mean(power)],'b*')

