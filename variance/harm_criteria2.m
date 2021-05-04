%%---------------------------------------------------------------------------------------------------
% function [No_Harm] = harm_criteria2(xt,xv,freq_sel,disp_error,ratio,freq_dat)
% Author: Maarten Buijsman
% Place: NIOZ, USM
% Date 1st version: 25-02-2004
% Date modified:    05-04-2021
% Description: 
% Determines minimum duration based on the Rayleigh's criterion,
% determines highest frequency based on the Nyquist criterium, and 
% determines lowest frequency
% input: vector time series xt [1]x[time] in decimal yearday, 
%        velocity series xv [1]x[time], vector freq_sel with number of
%        tidal constituents defined below, and disp_error (=1 YES display, =0 NO )
%        ratio ratio of data without NaNs to total data including NaNs
% output: No_Harm (=1 vector=bad, =0 vector=OK)
%
% use frequencies from frequencies_L2, (list of all frequencies)
% fq_str = string_freq(i,freq_dat), (produces strings with constituent name)

function [No_Harm] = harm_criteria2(xt,xv,freq_sel,disp_error,ratio,freq_dat)

%%-- frequencies [rad/day]
if strcmp(freq_dat,'L0')
    [fq] = frequencies;
elseif strcmp(freq_dat,'L1')
    [fq] = frequencies_L1;
elseif strcmp(freq_dat,'L2')
    [fq] = frequencies_L2;
end

%% find matching frequencies [rad/day]
fq2   = fq(freq_sel);	

%%-- check frequency resolution and determine min. duration time series for every cell
%%-- Reighly Criterium => minimum difference in frequency resolution: df = |f1-f2| = 1/T (T is duration time series)
fq_ch = fq2./2/pi; %%-- frequency [cycles/day]

if length(fq_ch) > 1

    [fq_ch_sort,index_sort] = sort(fq_ch);
    [min_diff,index_min] = min(diff(fq_ch_sort));
    min_T = 1/(min_diff);
 
    %% indices of neighbouring frequencies with smallest difference
    I1 = freq_sel(index_sort(index_min)); I2 = freq_sel(index_sort(index_min+1));
    %disp(['Bin ',num2str(i),'; Min. dur.: ',num2str(min_T),' days, due to fq: ',num2str(I1),' and ',num2str(I2)]);

    %%-- remove NANS
    le_xt = length(xt);
    index_NAN1  = find(isnan(xv)); index_NAN2  = find(isnan(xt));
    index_NAN = [index_NAN1 index_NAN2];
    xv(index_NAN) = [];
    xt(index_NAN) = [];

    fake_time = 0.50;
    if length(xt)/le_xt <= ratio | length(xt) == 1
        No_Harm = 1;        % i.e. less than 75 % of data available
        if disp_error == 1; disp(['WARNING - Ratio is ',int2str(length(xt)/le_xt*100),'%, less than ',int2str(ratio*100),'% of data available']); end
    else
        dt = diff(xt); %%-- difference between subsequent time steps
        if max(fq_ch) > 1/(2*min(dt)) %%-- Nyquist frequency (OK)
            min(dt)
            [dum,Imax] = max(fq_ch);
            Imax = freq_sel(Imax); %% locate correct one in freq_sel
            No_Harm = 1;
            if disp_error == 1; 
                disp(['Nyquist; freq. ',string_freq(Imax,freq_dat),'(',num2str(Imax),'); Min. dur.: ',num2str(1/max(fq_ch)),' days']);
            end
            %elseif fake_time > (xt(end) - xt(1)) %trick!! in the case of M2 M6 and M4
        elseif 1/min(fq_ch) > (xt(end) - xt(1)); %%-- lowest frequency (OK)
            [dum,Imin] = min(fq_ch);
            Imin = freq_sel(Imin); %% locate correct one in freq_sel
            No_Harm = 1;
            if disp_error == 1; 
                disp(['Lowest freq; freq. ',string_freq(Imin,freq_dat),'(',num2str(Imin),'); Min. dur.: ',num2str(1/min(fq_ch)),' days']);
            end
            %elseif length(fq2) > 1 & fake_time > (xt(end) - xt(1)); %trick!! in the case of M2 M6 and M4
        elseif length(fq2) > 1 & min_T > (xt(end) - xt(1));  %%-- Rayleigh, separation of frequencies, max duration (OK)
            No_Harm = 1;
            if disp_error == 1; 
                disp(['Freq resolution; freq. ',string_freq(I1,freq_dat),'(',num2str(I1),') and ',string_freq(I2,freq_dat),'(',num2str(I2),')']);
                disp(['     Min. dur.: ',num2str(min_T),' days']);
            end
        else
            No_Harm = 0;
        end
    end

elseif length(fq_ch) == 1

    le_xt = length(xt);
    index_NAN  = find(isnan(xv)); %%-- remove NANS
    xv(index_NAN) = [];
    xt(index_NAN) = [];

    if length(xt)/le_xt <= ratio
        No_Harm = 1;        % i.e. less than ratio*100% of data available
        if disp_error == 1; disp(['Ratio is ',int2str(length(xt)/le_xt*100),'%, less than ',int2str(ratio*100),'% of data available']); end
    else
        dt = diff(xt); %%-- difference between subsequent time steps
        if max(fq_ch) > 1/(2*min(dt)) %%-- Nyquist frequency (OK)
            No_Harm = 1;
            if disp_error == 1; disp(['Nyquist']); end
            %    elseif freq_sel == 28 & 0.516 > (xt(end) - xt(1)) % trick to bypass just too short time-series problem east/west 08/07/03
            %    elseif freq_sel == 28 & 0.49 > (xt(end) - xt(1)) % trick to bypass just too short time-series problem helsdeur stuff
        elseif 1/min(fq_ch) > (xt(end) - xt(1)); %%-- lowest frequency (OK)
            No_Harm = 1;
            if disp_error == 1; disp(['lowest freq']); end
        else
            No_Harm = 0;
        end
    end
end