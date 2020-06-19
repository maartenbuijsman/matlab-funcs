%% function [xp,yp] = draw_poly(closeflg,clr)
%  MCB, NRL, 2013-03-22
%  draws and prints polygon
%  stops with data input when "enter" is keyed
%  closeflg=1 closes polygon
%  if closeflg is omitted, polygon is closed (default)
%  default color of coordinates is white

function [xp,yp] = draw_poly(varargin)

% in default case polygon is always closed
if nargin == 0
    closeflg = 1;
    clr = 'w'
else
    closeflg = varargin{1};
    clr      = varargin{2};    
end


flg = 0
k=0;
while flg == 0
    % when enter is hit a and b are empty, and flg=1
    [a,b] = ginput(1);
    
    if ~isempty(a) % is empty when enter is hit
        k=k+1;
        
        holder
        plot(a,b,'ks','MarkerfaceColor',clr) % plot points
        
        if k>1; % connect with line
            plot([a c],[b d],['k--'])
        end
        c=a; d=b; % store points
       
        xp(k) = a;
        yp(k) = b;
    else
        flg = 1;
    end
end

% close it
if closeflg==1
    xp = [xp xp(1)];
    yp = [yp yp(1)];
    plot(xp(end-1:end),yp(end-1:end),['k--'])
end

% display values
disp(['xp = [' sprintf('%12.6f',xp) '];']);
disp(['yp = [' sprintf('%12.6f',yp) '];']);

return


