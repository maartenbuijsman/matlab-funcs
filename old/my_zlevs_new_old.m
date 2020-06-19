function [z,Cs] = my_zlevs_new(h,zeta,theta_s,theta_b,hc,N,type,scoord,alpha,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function z = my_zlevs_new(h,zeta,theta_s,theta_b,hc,N,type,scoord,alpha,beta)
%
%  this function compute the depth of rho or w points for ROMS
%
%  On Input:
%
%    type    'r': rho point 'w': w point
%    scoord     : 'old' (Song, 1994),
%                 'new' (Sasha, 2006)
%                 'bot' bottom stretching included (2008)
%    alpha,beta : used for 'bot'-type s-coordinate
%
%  On Output:
%
%    z       Depths (m) of RHO- or W-points (3D matrix).
%
%  Further Information:
%  http://www.brest.ird.fr/Roms_tools/
%
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven
%  e-mail:Pierrick.Penven@ird.fr
%
%  modified by Yusuke Uchiyama, UCLA, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 8
    error('Not enough input arguments')
elseif nargin < 9 & (scoord=='bot' | scoord=='BOT')
    disp('option bot is chosen, but not enough arguments. default values are taken.')
    alpha=0; beta=1;
end

%
[M, L] = size(h);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
if type=='w'
    sc = ((0:N) - N) / N;
    N = N + 1;
else
    sc=((1:N)-N-0.5) / N;
end

if (scoord=='bot' | scoord=='BOT');
    % new s-coordinate allowing smooth bottom refinement
    % alpha=-1: return to pure surface s-coord; -1 < alpha <\infty: bottom refinement
    % beta:
    if (theta_b>0);
        x = sc+1.0;
        wgt  = x.^alpha./beta.*(alpha+beta-alpha.*x.^beta);
        csrf = (1.0-cosh(theta_s*sc))./(cosh(theta_s)-1.0);
        cbot = sinh(theta_b*x)./sinh(theta_b)-1.0;
        Cs   = wgt.*csrf+(1.0-wgt).*cbot;
    else;
        Cs   = (1.0-cosh(theta_s*sc))./(cosh(theta_s)-1.0);
    end;
else;
    % for 'old' and 'new' s-coordinate
    cff1 = 1./sinh(theta_s);
    cff2 = 0.5/tanh(0.5*theta_s);
    Cs = (1.-theta_b) * cff1 * sinh(theta_s * sc)...
        + theta_b * (cff2 * tanh(theta_s * (sc + 0.5)) - 0.5);
end;
%
% Create S-coordinate system: based on model topography h(i,j),
% fast-time-averaged free-surface field and vertical coordinate
% transformation metrics compute evolving depths of of the three-
% dimensional model grid.
%
z=zeros(N,M,L);
if (scoord=='old' | scoord=='OLD')
    %disp('--- using old s-coord')
    hinv=1./h;
    cff=hc*(sc-Cs);
    cff1=Cs;
    cff2=sc+1;
    for k=1:N
        z0=cff(k)+cff1(k)*h;
        z(k,:,:)=z0+zeta.*(1.+z0.*hinv);
    end
else
    %disp('--- using new s-coord')
    hinv=1./(h+hc);
    cff=hc*sc;
    cff1=Cs;
    for k=1:N
        z(k,:,:)=zeta+(zeta+h).*(cff(k)+cff1(k)*h).*hinv;
    end
end

return

