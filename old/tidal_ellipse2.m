%%---------------------------------------------------------------------------------------------------
%% function [Amax,Amin,Ec,phi,psi,Rplus,Rmin,phiplus,phimin]=tidal_ellipse2(a_e,b_e,a_n,b_n);
%% Maarten Buijsman, NIOZ, 30-9-2003
%% Date 1st version: 16-7-2003 
%% Description: 
%% produces 9 parameters Amax,Amin,Ec,psi(inclination),phi(phase),Rplus,Rmin,phimin,phiplus
%% excentricity Ec, inclination psi, phase of maximum velocity phi,
%% and the maximum and minimum amplitudes Rmin and Rmax using coefficents a and b and 
%% uses Soulsby(1990) amd Prandle(1982)
%% Positive Amin and Ec indicate counterclockwise (cyclonal rotation)
%% note that phase phi is between t=0 and first max. velocity
%%---------------------------------------------------------------------------------------------------
function [Amax,Amin,Ec,phi,psi,Rplus,Rmin,phiplus,phimin]=tidal_ellipse2(a_e,b_e,a_n,b_n);

%% test
% a_e = gs2(1).a_e; b_e = gs2(1).b_e;
% a_n = gs2(1).a_n; b_n = gs2(1).b_n;

R1a = 1/2*(a_e+b_n); R1b = 1/2*(a_n-b_e);
R2a = 1/2*(a_e-b_n); R2b = 1/2*(a_n+b_e);

Rplus = sqrt(R1a.^2+R1b.^2);
Rmin  = sqrt(R2a.^2+R2b.^2);
phiplus = atan2(R1b,R1a);
phimin = atan2(R2b,R2a);

%% correct for -quadrants
%% check why this does not work, look at prandle
% If = find(phiplus<0); phiplus(If) = phiplus(If)+2*pi;
% If = find(phimin<0);  phimin(If) =  phimin(If)+2*pi;

Amax = Rplus+Rmin;
Amin = Rplus-Rmin;
phi = 1/2*(phimin-phiplus);
psi = 1/2*(phimin+phiplus);
Ec = (Rplus-Rmin)./(Rplus+Rmin);

% does this work??
phi(phi<0) = phi(phi<0)+pi;
psi(phi<0) = psi(phi<0)+pi;

phi(phi>2*pi) = phi(phi>2*pi)-2*pi;
psi(psi>2*pi) = psi(psi>2*pi)-2*pi;

