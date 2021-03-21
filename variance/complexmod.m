%% function [ACCW,ACW,PHICCW,PHICW,uf,vf] = complexmod(t,u,v,T,plotfig)
% Maarten Buijsman, USM, 2021-3-17
% Complex modulation to separate CW and CCW motions for prescribed frequency
%
% based on Emery and Thompson 3rd ed, page 557
% solve D*z=y (e.g. A*x=B)
%
% input: time series t, velocities u and v, prescribed period T
%        t,u,v should be 1D vectors 
%        t and T have the same time dimension (hours or days)
% output: Amplitudes ACCW and ACW, and phases PHICCW and PHICW in radians

function [ACCW,ACW,PHICCW,PHICW,uf,vf] = complexmod(t,u,v,T,plotfig);

% %% test -------------------------------------------
% t = 0:100; %hours
% T = 20;
% om = 2*pi/T;
% 
% % CW
% Au = 2;
% Av = 1;
% u = Au*cos(om*t);
% v = Av*cos(om*t-pi/2);
% figure; plot(u,v); axis equal
% 
% plotfig = 1
% % test -------------------------------------------

% check if row vectors
if size(u,1)>1; u=u'; end
if size(v,1)>1; v=v'; end
if size(t,1)>1; t=t'; end

om = 2*pi/T; %rad/hour

% create y matrix with stacked u and v
y = [u v]';

lu = length(u);

% create D matrix with 4 columns 
D=[];
D(:,1) = [ cos(om*t)  sin(om*t) ]';
D(:,2) = [-sin(om*t)  cos(om*t) ]';
D(:,3) = [ cos(om*t) -sin(om*t) ]';   
D(:,4) = [ sin(om*t)  cos(om*t) ]';  

%% solve it
z = D\y;

% map it to E&T variables
ACP = z(1);
ASP = z(2);
ACM = z(3);
ASM = z(4);

% create the fit yf
yf = D*z;

uf = yf(1:lu);
vf = yf(lu+1:end);
    
%% extract Amplitudes ACCW and ACW, and phases PHICCW and PHICW in radians
ACCW = sqrt(ASP^2+ACP^2);
ACW  = sqrt(ASM^2+ACM^2);

PHICCW = atan2(ASP,ACP);
PHICW  = atan2(ASM,ACM);

%% plot output if asked for
if plotfig
    figure
    subplot(2,1,1)
    plot(t,u,'r-',t,uf,'b--')

    subplot(2,1,2)
    plot(t,v,'r-',t,vf,'b--')
end

