function d=bandpass_new(c,T1,T2,delt,N)
%
% d=bandpass_new(c,T1,T2,delt,N)
%
% MCB, USM, 2021-3-19
%
% bandpass a time series with a 2*Nth order butterworth filter
%
% c = input time series
% T1 = smallest period of filter (left)
% T2 = largest period of filter (right)
% delt = sampling interval of data
% N    = yields 2*Nth order butterworth filter 

% convert periods to cutoff frequencies
fcutlow  = 1/T2;
fcuthigh = 1/T1;

fnq=1/(2*delt);             % Nyquist frequency
Wn=[fcutlow fcuthigh]/fnq;  % butterworth bandpass non-dimensional frequency
[b,a]=butter(N,Wn);         % construct the filter
d=filtfilt(b,a,c);          % zero phase filter the data

return;

