%%---------------------------------------------------------------------------------------------------
%%-- Program: figA4L.m
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 01-10-2003
%%-- Date 1st programmed: 01-10-2003
%%-- Description: 
%%-- makes landscape A4 size figure
%%---------------------------------------------------------------------------------------------------
function f1 = figA4L;

%f1 = figure('paperunits','centimeters','PaperPosition',[0.5 0.5 28.6774 19.984],'Position',[296 50 800 600],'PaperType','A4','PaperOrientation','landscape');
%f1 = figure('PaperUnits','centimeters','PaperPosition',[0.5 0.5 28.6774 17],'Position',[296 100 1200 700],'PaperType','A4','PaperOrientation','portrait');
%f1 = figure('PaperUnits','normalized','PaperPosition',[0.02 0.02 0.96 0.96],'Position',[296 50 950 800],'PaperType','usletter','PaperOrientation','landscape');
f1 = figure('PaperUnits','inches','PaperPosition',[0.22 0.22 10.5600 8.0600],'Position',[100 350 1000 600],'PaperSize',[11 8.5]);

