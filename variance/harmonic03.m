%% function [Ao,a,b,fq2,coef_det,sigma,fit]=harmonic03(xt,xv,freq_sel,cel_num,plot_fig,freq_dat)
%  Author: Maarten Buijsman
%  Place: NIOZ, USM
%  Date modified:    04-05-2021
%  Date 1st version: 21-03-2003
%
%  Least squares harmonic analysis for timeseries [depth vs time]
%  input: xt vector times [1]*[nt] in decimal yearday, 
%         xv matrix velocity time series [n]*[nt], 
%         freq_sel vector with number of tidal constituents defined as in frequencies_L2.m 
%         cel_num index of depth level of xv; if only vector is used cel_num=1
%         plot_fig = 1 if figure of time series and fit needs to be displayed
%         freq_dat = 'L2' if frequencies_L2 is used to obtain frequencies
%
%  output: the time-mean Ao, constituent coefficients a and b 
%          as in Ao + sum(n){a(n)*cos[omega(n)*t]+b(n)*sin[omega(n)*t]} 
%          coef_det, the coefficient of determination R^2 = 1-diff^2/tot^2
%          sigma, the standard deviation of the residual
%          fit, fitted time series with the same size as xv 
%  
%  test from Emery and Thomson, p401, fitting to k1 and m2 freq_sel = [21 46]
%  xt=[1:1:32]./24 %  time in days
%  xv=[1.97 1.46 0.98 0.73 0.67 0.82 1.15 1.58 2.0 2.33 2.48 2.43 2.25 2.02 1.82 1.72 1.75 1.91 2.22 2.54 2.87 3.1 3.15 2.94 2.57 2.06 1.56 1.13 0.84 0.73 0.79 1.07]

function [Ao,a,b,fq2,coef_det,sigma,fit]=harmonic03(xt,xv,freq_sel,cel_num,plot_fig,freq_dat)

%% frequencies [rad/day]
if strcmp(freq_dat,'L0')
    [fq] = frequencies;
elseif strcmp(freq_dat,'L1')
    [fq] = frequencies_L1;
elseif strcmp(freq_dat,'L2')
    [fq] = frequencies_L2;
elseif strcmp(freq_dat,'L2test')
    [fq] = frequencies_L2test;
end

%% find matching frequencies
fq2   = fq(freq_sel);	%[rad/day]

%% number of variables in signal
%  z(t) = mean + acos(fq*t) + bsin(fq*t)
%  Per component 2 variables a and b, plus 1 variable for the mean current
ntide = 2*length(freq_sel)+1;		

%% make sure it is a row matrix
%if size(xv,1)>size(xv,2); 
%    xv = xv'; 
%end

%%-- number of depth cells
nbins = size(xv,1);	

if plot_fig == 1; figure; end

%% method used is the least squares analysis presented in numerical recepies chapter 14.3
%  declare matrices to make it FASTER!
x2  = zeros(ntide,nbins); %%-- definine matrix [number of variables * depth bins]
Ao  = zeros(nbins,1)*NaN;                coef_det=Ao; sigma=Ao;
a   = zeros(length(freq_sel),nbins)*NaN; b=a;
fit = zeros(nbins,length(xt))*NaN;       data = fit;
for x4 = 1:nbins;	     %%-- Loop over all bins (Harmonic analysis for all depths)
    xv2 = xv(x4,:)';	 %%-- current velocities of bin number x4 (xv2 = cell(nbins) x time)
    xt2 = xt;		     %%-- time vector
    x5  = find(~isnan(xv2));	%%-- remove NaN's verwijderen; find indices of elements without NaN for all time steps for one cel(x4)
    xv2 = xv2(x5);       %%-- calc. matrices for velocity and time with new indices (without NaNs)
    xt2 = xt2(x5);
   
    if length(xv2) > ntide;	%% Check if number of data points is smaller than number of a's and b's (too few => not a good fit possible)
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
        %  see Emery and Thomson, p.235
        xv_harm      = har2*x2(:,x4);
        SSR          = (xv_harm-mean(xv_harm))'*(xv_harm-mean(xv_harm));        
        SST          = (xv2-mean(xv2))'*(xv2-mean(xv2));       
        SSE          = (xv2-mean(xv2)-(xv_harm-mean(xv_harm)))'*(xv2-mean(xv2)-(xv_harm-mean(xv_harm)));
        %coef_det(x4) = SSR/SST;
        coef_det(x4) = 1-SSE/SST;   %this one is better according to Wikipedia!
        sigma(x4)    = sqrt(SSE/(x6-ntide));
        
        %% plot for value x4 (in case bottom-up, 1st cell above the bottom is 2nd cell!)
        if plot_fig == 1; 
           if x4 == cel_num; 
                plot(xt2,xv2,'b.-.'); hold; %data
                plot(xt2,xv_harm,'r.-');    %fit
           end
        end
        
        %% fit will be exported
        fit(x4,x5)  = xv_harm';
        data(x4,x5) = xv2';
        
    else  %% in case length(xv2) smaller than ntide no harm.ana.
        %disp(['Bin ',num2str(x4+1),' skipped'])
        Ao(x4)              = NaN;  %% residual current   
        a(1:(ntide-1)/2,x4) = NaN;	%% a-value in acos(fq*t)
        b(1:(ntide-1)/2,x4) = NaN;	%% b-value in bsin(fq*t)
        coef_det(x4)        = NaN;
        sigma(x4)           = NaN;        
    end;
    clear x5 xt2 xv2 ;
end;




















