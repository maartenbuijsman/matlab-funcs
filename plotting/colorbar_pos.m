%% function h1 = colorbar_pos(fighandle,figpos,figbarspace,wbar,hbar,vertoffbar)
%% add colorbar to figure
%% without rescaling the figure
%% hbar and vertoffbar = fraction of figure height
%% figbarspace, wbar are relative to paper
%% 
%% example
%% figbarspace=0.025; wbar=0.025; hbar=1; vertoffbar=0;
%% h1 = colorbar_pos(fh,pos(1,:),figbarspace,wbar,hbar,vertoffbar);


function h1 = colorbar_pos(fighandle,figpos,figbarspace,wbar,hbar,vertoffbar)

h1=colorbar('position',[figpos(1)+figpos(3)+figbarspace figpos(2)+vertoffbar*figpos(4) wbar hbar*figpos(4)]);
set(fighandle,'position',figpos);
