%%---------------------------------------------------------------------------------------------------
%% function: tidal_ellipse.m
%% Author: Maarten Buijsman
%% Place: NIOZ
%% Date: 30-9-2003
%% Date 1st version: 16-7-2003 
%% Function links: 
%% Description: 
%% program calculates the excentricity Ec, inclination psi, phase of maximum velocity phi,
%% and the maximum and minimum amplitudes Rmin and Rmax using coefficents a and b and 
%% angular frequency omega from the harmonic analysis 
%% 30-09-03 included correction for negative (clockwise) or positive (anti-clockwise) eccentricity
%%---------------------------------------------------------------------------------------------------
function [Amax,Amin,Ec,phi,psi]=tidal_ellipse(a_e,b_e,a_n,b_n,omega)
%% calculate tidal ellips according to Prandle, 1982

%% circel R1 
R1_e = 1/2*(a_e + b_n); R1_n = 1/2*(a_n - b_e);
R1_abs = sqrt(R1_e.^2+R1_n.^2);
g1 = atan2(R1_n,R1_e); g1_deg = g1 * 180/pi;
%% circel R2
R2_e = 1/2*(a_e - b_n); R2_n = 1/2*(a_n + b_e);
R2_abs = sqrt(R2_e.^2+R2_n.^2);
g2 = atan2(R2_n,R2_e); g2_deg = g2 * 180/pi;

mu = R1_abs./R2_abs;
Ec = (mu-1)./(mu+1);
Amax = R1_abs + R2_abs;
Amin = abs(R1_abs - R2_abs).*sign(Ec);

tmax = (g2-g1)/(2*omega);
phi = omega.*tmax;   phi_deg = phi .* 180/pi; % phase of max velocity relative to t0
psi = 1/2*(g2+g1);   psi_deg = psi .* 180/pi;  % inclination relative to the x-axis
