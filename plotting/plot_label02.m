%%---------------------------------------------------------------------------------------------------
%%-- function plot_label02(FontSize,LW,title_lab,x_lab,y_lab);
%%-- Author: Maarten Buijsman
%%-- Place: NIOZ
%%-- Date: 10-1-2003
%%-- Description: 
%%-- plots axis labels 
%%-- input: [FontSize,title_lab,x_lab,y_lab]
%%-- output: - 
%%---------------------------------------------------------------------------------------------------
function plot_label02(FontSize,LW,title_lab,x_lab,y_lab);

h1 = title(title_lab);  set(h1,'FontSize',FontSize,'FontWeight','normal');
h2 = xlabel(x_lab);     set(h2,'FontSize',FontSize);
h3 = ylabel(y_lab);     set(h3,'FontSize',FontSize);
set(gca,'FontSize',FontSize);
set(gca,'Layer','top');
set(gca,'LineWidth',LW);
