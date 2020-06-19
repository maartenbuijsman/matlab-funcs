function [puv,quv,cw,ccw] = ff_spec_rot(u, v);
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

% individual components fourier series
% units: u_unit
fu = fft(u); fv = fft(v);

% autospectra of the scalar components
% units: u_unit^2
pu = fu.*conj(fu); pv = fv.*conj(fv);

% cross spectra
% units: u_unit^2
puv =  real(fu).*real(fv) + imag(fu).*imag(fv);

% quadrature spectra
% units: u_unit^2
quv = -real(fu).*imag(fv) + real(fv).*imag(fu);

% rotatory components
% units: u_unit^2
cw   = (pu+pv-2*quv)./8;
ccw  = (pu+pv+2*quv)./8;
