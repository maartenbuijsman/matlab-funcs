function [ Pxx, Pyy, Pxy, coh, pha, freq ] = func_coherence_win( f1, f2, win, n_overlap, nfft, Fs )

%
% Computing power spectrum, cross spectrum, coherence and phase
%
% All input parameters are equivalent to csd or the other spectrum 
% related function in MATLAB.
%
% This program use csd.m in Signal Processing Toolbox.
%
% Input
%	f1 and f2	input data
%	nfft		number of data for FFT; default = 256
%	Fs          sampling frequency; default = 1 Hz
%	win		(=window) filter vector, hanning(nfft/2)
%	n_overlap	number of overlap for smoothing
%
% Output
%	Pxx		spectrum of f1
%	Pyy		spectrum of f2
%	Pxy		cross spectrum between f1 and f2
%	coh		coherence
%	pha		phase
%	freq		frequency vector
%
% Date: 2005-08-02
% From: Vitaliy Marchenko (marchenk@mail.eecis.udel.edu) 
% Rating:  
% Comments: check it with two aboslutely identical sinusoidal waves of any frequenciy (let say) - 50 hz ... - and you will be surprized... - it doesn't work as Matlab standart function - 'cohere' or 'mscohere'(for version 7.0)  
%  
%  Date: 2005-04-16
% From: Philip Orton (orton@ldeo.columbia.edu) 
% Rating:  
% Comments: I agree. Simple, but just what I was looking for, since I wasn't aware of how to compute phase. A minor improvement would be to allow the user to ignore the filt and n_overlap inputs, just as one can do with the csd.m function.  
%  
%  Date: 2004-04-02
% From: A Zaki (zaki@egr.uri.edu) 
% Rating:  
% Comments: Good program..  
%
% Date: 22 Jun 2006
% From: Jon Fram   	
%I got it to produce the same output as mscohere by replacing csd with cpsd.
%I also wasn't clear how to calculate phase. func_coherence.m makes that clear to me. Thanks.

% default
[Pxx,freq] = cpsd( f1, f1, win, n_overlap, nfft, Fs);
[Pyy,freq] = cpsd( f2, f2, win, n_overlap, nfft, Fs);
[Pxy,freq] = cpsd( f1, f2, win, n_overlap, nfft, Fs);

Kxy  = real( Pxy );
Qxy  = imag( Pxy );
coh  = Pxy.*conj(Pxy)./(Pxx.*Pyy);
pha  = atan2( Qxy, Kxy );




