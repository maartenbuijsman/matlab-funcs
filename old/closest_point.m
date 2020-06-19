%% Ipnt = closest_point(x,xpnt)
%  find grid point x closest to xpnt
%  MCB, NRL, 2013-01-16

function Ipnt = closest_point(x,xpnt)

[d,Ipnt] = min(abs(x-xpnt));

