%% function [uu]=plot_sinus(cc,plotflg)
%% plots sinusoid given complex vector cc

function uu=plot_sinus(cc,plotflg);

clr = {'k' 'b' 'r' 'y' 'g' 'c' 'm'};
tt=linspace(0,2*pi,50);

for i=1:length(cc)
    uu(i,:) = abs(cc(i))*cos(tt-angle(cc(i)));

    if plotflg 
        if i==1; figure; end
        
        plot(tt,uu(i,:),char(clr(i)))
        holder
    end
    
end

axis tight
title(char(clr(1:i))')

 
