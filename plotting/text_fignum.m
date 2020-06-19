%% function handle = text_fignum(str,offhorz,offvert,horloc,vertloc,fs,bckclr)
%% MCB, UCLA, 2009-05-01
%% prints figure label on left top corner of figure
%% offhorz: is relative position relative to width of x axis; -0.1 is good value
%% offvert: is relative position relative to width of y axis; 0 is good value
%% vertloc: is location on 'left' or 'right' of plot
%% horloc:  is location on 'top' or 'bottom' of plot
%% fs:      is fontsize
%% bckclr:  is background color (t=transclucent) 

function handle = text_fignum(str,offhorz,offvert,horloc,vertloc,fs,bckclr);

xlim=get(gca,'xlim'); dx=diff(xlim);
ylim=get(gca,'ylim'); dy=diff(ylim);

%offr = -0.1; 
offh = offhorz*dx;
offv = offvert*dy;

if strcmp(vertloc,'top')
    if strcmp(horloc,'left')
        handle = text(xlim(1)+offh,ylim(2)+offv,['\bf' str '\rm'],'VerticalAlignment','top','HorizontalAlignment','Le','fontsize',fs);
    elseif strcmp(horloc,'right')
        handle = text(xlim(2)+offh,ylim(2)+offv,['\bf' str '\rm'],'VerticalAlignment','top','HorizontalAlignment','Ri','fontsize',fs);
    end
elseif strcmp(vertloc,'bottom')
    if strcmp(horloc,'left')
        handle = text(xlim(1)+offh,ylim(1)+offv,['\bf' str '\rm'],'VerticalAlignment','bottom','HorizontalAlignment','Le','fontsize',fs);
    elseif strcmp(horloc,'right')
        handle = text(xlim(2)+offh,ylim(1)+offv,['\bf' str '\rm'],'VerticalAlignment','bottom','HorizontalAlignment','Ri','fontsize',fs);
    end
end

if ~strcmp(bckclr,'t')
    set(handle,'BackgroundColor',bckclr)
end

return