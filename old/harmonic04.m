%%---------------------------------------------------------------------------------------------------
%%-- function: [Ao,a,b,fq2,coef_det,sigma,data,fit] = harmonic04(xt,xv,freq_sel,cel_num,plot_fig,freq_dat)
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date modified:    24-12-2003
%%-- Date 1st version: 21-3-2003
%%-- Description: 
%%-- Harmonic analysis for velocity vectors
%%-- applies least squares technique
%%-- input: vector time series xt [1]x[time] in decimal yearday, velocity series xv [n]x[time], 
%%-- and vector freq_sel with number of tidal constituents defined below:
%%-- output: the residual Ao, coefficients a and b from Ao + sum(n){a(n)*cos[omega(n)*t]+b(n)*sin[omega(n)*t]} 
%%-- and the coefficient of determination (r^2) and the standard deviation
%%-- frequencies in [rad/day]:
%%-- Sa(1);    Ssa(2);  Mm(3);    Mf(4);   2Q1(5);   sig1(6);  q1(7);    ro1(8);   O1(9);   tau1(10); M1(11);   NO1(12);
%%-- x1(13);   pi1(14); P1(15);   S1(16);  K1(17);   psi1(18); phi1(19); the1(20); J1(21)   OO1(22);  eps2(23); 2N2(24); 
%%-- mu2(25);  N2(26);  nu2(27);  M2(28);  lab2(29); L2(30);   T2(31);   S2(32);   R2(33);  K2(34);   ksi2(35); eta2(36); 
%%-- M3(37);   SO1(38); OQ2(39);  OP2(40); MKS2(41); 2SM2(42); MO3(43);  SO3(44);  MK3(45); SK3(46);  MN4(47);  M4(48);
%%-- SN4(49)); MS4(50); MK4(51);  S4(52);  SK4(53);  2MN6(54); M6(55);   MSN6(56); 2MS6(57);SMK6(58); 2SM6(59); MSK6(60));
%%-- 3MN8(61); M8(62);  2MSN8(63);3MS8(64);2MS8(65); 2MSK8(66); M2S2 beat(68) 
%%-- 
%%-- test from Emery and Thomson, p401, fitting to k1 and m2 freq_sel = [17 28]
%%-- xt=[1:1:32]./24 %%-- time in days
%%-- xv=[1.97 1.46 0.98 0.73 0.67 0.82 1.15 1.58 2.0 2.33 2.48 2.43 2.25 2.02 1.82 1.72 1.75 1.91 2.22 2.54 2.87 3.1 3.15 2.94 2.57 2.06 1.56 1.13 0.84 0.73 0.79 1.07]
%%-- 30-09-03 added position cel_num for plotting
%%-- 24-12-03 returns the data and the fit matrix
%%-- 09-05-05 choice for plotting
%%-- 12-05-29 avoids zero time series => in that case output is zero
%%---------------------------------------------------------------------------------------------------
function [Ao,a,b,fq2,coef_det,sigma,data,fit]=harmonic04(xt,xv,freq_sel,cel_num,plot_fig,freq_dat)

%%-- frequencies [rad/day]
if strcmp(freq_dat,'L0')
    [fq] = frequencies;
elseif strcmp(freq_dat,'L1')
    [fq] = frequencies_L1;
elseif strcmp(freq_dat,'L2')
    [fq] = frequencies_L2;
end

%% find matching frequencies
fq2   = fq(freq_sel);	%[rad/day]

%% number of variables in signal
%% z(t) = mean + acos(fq*t) + bsin(fq*t)
%% Per component 2 variables a and b, plus 1 variable for the mean current
ntide = 2*length(freq_sel)+1;		

%% make sure it is a row matrix
%if size(xv,1)>size(xv,2); 
%    xv = xv'; 
%end

%%-- number of depth cells
nbins = size(xv,1);	

if plot_fig == 1; figure; end

%%-- method used is the least squares analysis presented in numerical recepies chapter 14.3
%% declare matrices to make it FASTER!
x2  = zeros(ntide,nbins); %%-- definine matrix [number of variables * depth bins]
Ao  = zeros(nbins,1)*0;                  coef_det=Ao; sigma=Ao;
a   = zeros(length(freq_sel),nbins)*0;   b=a;
fit = zeros(nbins,length(xt))*0;         data = fit;
for x4 = 1:nbins;	     %%-- Loop over all bins (Harmonic analysis for all depths)
    xv2 = xv(x4,:)';	 %%-- current velocities of bin number x4 (xv2 = cell(nbins) x time)
    xt2 = xt;		     %%-- time vector
    x5  = find(~isnan(xv2));	%%-- remove NaN's verwijderen; find indices of elements without NaN for all time steps for one cel(x4)
    xv2 = xv2(x5);       %%-- calc. matrices for velocity and time with new indices (without NaNs)
    xt2 = xt2(x5);
   
    if length(xv2) > ntide & std(xv2)~=0;	% Check if number of data points is smaller than number of a's and b's (too few => not a good fit possible)
                                            % and xv2 is not zero                
        x6    = length(xt2);	%% determine length of new variables; length of xt2 => time matrix
        har2  = zeros(x6,ntide); %% matrix A (NR 14.3) containing the eigenfunctions [x6 of rows and ntide of columns] 
        har2(:,1) = 1;		%% first column: residual current velocities => eigenfunction is 1
                          
        for x7  = 1:(ntide-1)/2; %% fill A: loop over all tidal components; cos and sin alternate columns; comp: 112233..; NR, p.510, eq.14.3.4 
            har2(:,2*x7)   = cos(fq2(x7)' .* xt2'); %% Eigenfunction cos(fq*t) for calc of a
            har2(:,1+2*x7) = sin(fq2(x7)' .* xt2'); %% Eigenfunction sin(fq*t) for calc of b
        end;
        [U,S,V]  = svd(har2,0);	%% Singular Value Decomposition (type help svd in matlab)
        for x8 = 1:ntide;
            if S(x8,x8) <= eps*x6; %% eps: machine-accuracy (can neglect this if loop; see pg 55 Num. Rec.)
                invw = 0;
            else invw = 1/S(x8,x8); %% use values of diagonal of S
            end;
            x2(:,x4) = x2(:,x4) + invw* (U(:,x8)'*xv2)*V(:,x8); %%-- All values of a en b are stored for every depth bin x4, summation over x2; NR, p.516, eq. 14.3.17 
        end;

        Ao(x4)  = x2(1,x4);	        %% residual current    
        a(:,x4) = x2(2:2:end,x4);	%% a-value in acos(fq*t)
        b(:,x4) = x2(3:2:end,x4);	%% b-value in bsin(fq*t)
           
        %% calculate coefficient of determination and standard deviation
        %% see Emery and Thomson, p.235
        xv_harm      = har2*x2(:,x4);
        SSR          = (xv_harm-mean(xv2))'*(xv_harm-mean(xv2));        
        SST          = (xv2-mean(xv2))'*(xv2-mean(xv2));       
        coef_det(x4) = SSR/SST;
        SSE          = (xv2-xv_harm)'*(xv2-xv_harm);
        sigma(x4)    = sqrt(SSE/(x6-ntide));
        
        %% plot for value x4 (in case bottom-up, 1st cell above the bottom is 2nd cell!)
        if plot_fig == 1; 
           if x4 == cel_num; 
                plot(xt2,xv2,'b+--'); hold; %data
                plot(xt2,xv_harm,'r.-');    %fit
           end
        end
        
        %% fit will be exported
        fit(x4,x5)  = xv_harm';
        data(x4,x5) = xv2';
        
    else  %% in case length(xv2) smaller than ntide no harm.ana.
        %disp(['Bin ',num2str(x4+1),' skipped'])
        Ao(x4)              = 0;  %% residual current   
        a(1:(ntide-1)/2,x4) = 0;	%% a-value in acos(fq*t)
        b(1:(ntide-1)/2,x4) = 0;	%% b-value in bsin(fq*t)
        coef_det(x4)        = 0;
        sigma(x4)           = 0;        
    end;
    clear x5 xt2 xv2 ;
end;




















