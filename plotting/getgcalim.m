%% function [xlim ylim] = getgcalim;
%% MVB, UCLA, Maarten Buijsman, 2010-02-09
%% get axis limits of current axis

function [xlim ylim] = getgcalim;

xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

disp(sprintf('%8.3f ',[xlim ylim]));
