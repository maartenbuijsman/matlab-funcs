%%---------------------------------------------------------------------------------------------------
%%-- function [u2,v2,theta_main] = eofr(u1,v1);
%%-- Author: Maarten Buijsman
%%-- Place: USM
%%-- Date: 09-12-2017
%%-- Description: 
%%-- Principal component analysis
%%-- Calculates primairy axis in radians (theta_main_rad) of the u1 and v1 data
%%-- Exports u2 and v2, velocities relative to the new axis
%%-- see Preisendorfer (1988) & Emery and Thomson (2001), p. 327
%%---------------------------------------------------------------------------------------------------

function [u2,v2,theta_main_rad] = eofr(u1,v1);

u = u1; %u1 and v1 still have NaNs
v = v1;

ff =  find(isnan(u)); %remove NaN
u(ff) = [];
v(ff) = [];
% gg = find(abs(u)>2.5);
% u(gg) = [];
% v(gg) = [];
% hh = find(abs(v)>2.5);
% u(hh) = [];
% v(hh) = [];

if length(u) > 1
    u_ac = u - mean(u); %remove the mean
    v_ac = v - mean(v); 
    
    Suu = 1/(length(u)-1)*sum(u_ac.^2);
    Svv = 1/(length(v)-1)*sum(v_ac.^2);
    Suv = 1/(length(u)-1)*sum(u_ac.*v_ac);

    numer = 2*Suv;
    denom = Suu-Svv;

    %%-- 9 cases
    if numer < 0 & denom < 0
        theta_main_rad = 1/2 * atan(numer/denom) - pi/2;
    elseif numer < 0 & denom == 0
        theta_main_rad = - pi/4;
    elseif numer < 0 & denom > 0
        theta_main_rad = 1/2 * atan(numer/denom);
    elseif numer == 0 & denom < 0
        theta_main_rad = pi/2;
    elseif numer == 0 & denom == 0
        theta_main_rad = 0;
    elseif numer == 0 & denom > 0
        theta_main_rad = 0;
    elseif numer > 0 & denom < 0
        theta_main_rad = 1/2 * atan(numer/denom) + pi/2;
    elseif numer > 0 & denom == 0
        theta_main_rad = pi/4;
    elseif numer > 0 & denom > 0
        theta_main_rad = 1/2 * atan(numer/denom);
    end

    theta_main = theta_main_rad*180/pi; %radians => degrees

    %%-- calculate new velocities relative to the primary and secondairy axis
    u2 = u1.*cos(theta_main_rad) + v1.*sin(theta_main_rad);
    v2 = -u1.*sin(theta_main_rad) + v1.*cos(theta_main_rad);
else    
    u2(1:length(u1)) = NaN;
    v2(1:length(u1)) = NaN;
    theta_main = NaN;
end



%for i=1:360
%    theta = i*180/pi;
%    Ssq(i) = Suu*cos(theta)*cos(theta)+2*Suv*sin(theta)*cos(theta)+Svv*sin(theta)*sin(theta);
%end
%plot(1:1:360,Ssq)
%max(Ssq)
%min(Ssq)
%theta = theta_main_rad;
%S11 = Suu*cos(theta)*cos(theta)+2*Suv*sin(theta)*cos(theta)+Svv*sin(theta)*sin(theta)
%theta = theta_main_rad + pi/2;
%S22 = Suu*cos(theta)*cos(theta)+2*Suv*sin(theta)*cos(theta)+Svv*sin(theta)*sin(theta)
%
%disp(['teller = ',num2str(numer)])
%disp(['noemer = ',num2str(denom)])    
%disp(['ratio = ',num2str(numer/denom)])    

%plot(u,v,'b.');
%axis equal
%grid
%hold off
%pause

%plot(u2,v2,'b.');
%hold on

%U(1) = 1;
%U(2) = -1;
%V(1) = tan(theta_main_rad)*U(1);
%V(2) = -V(1);
%h2 = quiver([0 0],[0 0],U,V);
%set(h2,'Color','r'); 
%grid
%hold off
%axis equal
%drawnow
%pause

%   clear all
%    i=1;
%    u=[];
%    v=[];
%    u(1:1000) = 1;
%    u(1001:2000) = 2;
%    u(2001:3000) = -1;
%    v(1:1000) = 2;
%    v(1001:2000) = 3;
%    v(2001:3000) =-2; 
