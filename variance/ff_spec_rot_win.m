function [puv,quv,cw,ccw,period,freq] = ff_spec_rot_win(t1,u1,v1,numwin);
% ff_spec_rot_win, modified by mcb, USM, 2020-7-1
%   t1 is time steps
%   numwin is number of windows (segments) that time series will be split in 
%       spectra will be averaged over numwin
%   period and freq are perod and frequency
%   u1 and v1 are x and y velocity components
% units of puv,quv,cw,ccw: u_unit^2 
%    in order to get  u_unit^2/Hz = u_unit^2 * dt
% ff_spec_rot.m -> compute the rotary spectra from u,v velocity components
%
% use:  [puv,quv,cw,ccw] = ff_spec_rot(u,v);
% input:
%        u-component of velocity
%        v-component of velocity
%
% output:
%        cw  - clockwise spectrum
%        ccw - counter-clockwise spectrum
%        puv - cross spectra
%        quv - quadrature spectra
%
% example:
%    [puv,quv,cw,ccw] = ff_spec_rot(u,v)
%
% other m-files required: fft.m
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     19-Jan-2002
% modified: 06-May-2009 (translated to english)
% https://code.google.com/p/ocefpaf-matlab/source/browse/ocfis/?name=38&r=ae9f51452c4fb49db8bb295a7cbf83be707006e0
%
% obs: based on J. Gonella Deep Sea Res., 833-846, 1972
%      definition: The spectral energy at some frequency can be
%      decomposed into two circulaly polarized constituents, one
%      rotating clockwise and other anti-clockwise
%

% %test
% clear all
% numwin = 3;
% t1=1:0.1:1000;
% u1 = sin(2*pi/12*t1);
% v1 = 2*cos(2*pi/12*t1-5);

% number of indices per numwin+1 segments
nt = length(t1);
inw = floor(nt/(numwin+1)); 

% get indices
is=1;
for i=1:numwin
    p(i).ii = is:2*inw+is-1;
    is = i*inw+1;
end

for i=1:numwin
    u = u1(p(i).ii);
    v = v1(p(i).ii);
    
    t = t1(p(i).ii);  
    
    % remove the mean!
    u = u - mean(u);
    v = v - mean(v);    

    % individual components fourier series
    fu = fft(u); fv = fft(v);
    
    % The first component is simply the sum of the data, and can be removed.
    fu(1)=[];
    fv(1)=[];
    
    % time releated ......
    n=length(fu); 
    dt = mean(diff(t));
    nyquist = 1/(2*dt);  %% highest frequency that can be resolved: 1/(2(dt)); 
    p(i).freq = (1:n/2)/(n/2)*nyquist;
    p(i).period=1./p(i).freq;
    
    % autospectra of the scalar components
    pu = fu.*conj(fu); pv = fv.*conj(fv);

    % cross spectra
    p(i).puv =  real(fu).*real(fv) + imag(fu).*imag(fv);

    % quadrature spectra
    p(i).quv = -real(fu).*imag(fv) + real(fv).*imag(fu);

    % rotatory components
    p(i).cw   = (pu+pv-2*p(i).quv)./8;
    p(i).ccw  = (pu+pv+2*p(i).quv)./8;

end

% finally do some averaging
puv=0; quv=0; cw=0; ccw=0;
for i=1:numwin
    puv = puv + p(i).puv(1:floor(n/2));
    quv = quv + p(i).quv(1:floor(n/2));
    cw  = cw  + p(i).cw(1:floor(n/2));
    ccw = ccw + p(i).ccw(1:floor(n/2));    
end
puv = puv/numwin;
quv = quv/numwin;
cw  = cw/numwin;
ccw = ccw/numwin;    

period=p(1).period;
freq=p(1).freq;

% figure
% plot(period,log10(cw))
% hold
% plot(period,log10(ccw),'k.-')

