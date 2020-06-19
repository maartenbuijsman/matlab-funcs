%% function [phi2,psi2] = straighten_phi_psi(phi2,psi2,topdown,ang);
%% Maarten Buijsman, UCLA, 2009-12-09
%% finds coinciding jumps in phi2 and psi2 from tidal_ellipse2.m
%% and removes than in both phi2 and psi2
%% INPUT:  vectors phi2 and psi2; topdown='surf', when doing "unwrap"ping from surface
%% and angle ang; when [] default ang = pi*3/4 is used 
%% OUTPUT:  vectors phi2 and psi2

function [phi2,psi2] = straighten_phi_psi(phi2,psi2,topdown,ang);

% %% test
% topdown='surf'; ang=[];
% phi2 = phi(:,Ip);
% psi2 = psi(:,Ip);
% 
% figure
% plot(phi2*180/pi,'k.-')
% hold;
% plot(psi2*180/pi,'r.-')
% %% test

%% default angle
%if isempty(ang); ang = pi*3/4; end
if isempty(ang); ang = pi/2*5/4; end

if strcmp(topdown,'surf')
    phi2 = phi2(end:-1:1);
    psi2 = psi2(end:-1:1);    
end

%% find jumps
Ij1 = find(abs(diff(phi2))>ang & abs(diff(psi2))>ang);

for k=1:length(Ij1)
    DP = phi2(Ij1(k)+1)-phi2(Ij1(k));
    phi2(Ij1(k)+1:end) = phi2(Ij1(k)+1:end) - pi*round(DP/pi);
    psi2(Ij1(k)+1:end) = psi2(Ij1(k)+1:end) + pi*round(DP/pi);
end

%% remove any jumps
phi2 = unwrap(phi2);
psi2 = unwrap(psi2);

%% reverse to bottup-up
if strcmp(topdown,'surf')
    phi2 = phi2(end:-1:1);
    psi2 = psi2(end:-1:1);    
end
