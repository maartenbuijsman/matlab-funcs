function [zi] = func_fast_griddata(x,y,z,xi,yi,method,tri,td)
%GRIDDATA Data gridding and surface fitting.
%   ZI = GRIDDATA(X,Y,Z,XI,YI) fits a surface of the form Z = F(X,Y) to the
%   data in the (usually) nonuniformly-spaced vectors (X,Y,Z). GRIDDATA
%   interpolates this surface at the points specified by (XI,YI) to produce
%   ZI.  The surface always goes through the data points. XI and YI are
%   usually a uniform grid (as produced by MESHGRID) and is where GRIDDATA
%   gets its name.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it specifies
%   a matrix with constant rows.
%
%   [XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) also returns the XI and YI formed
%   this way (the results of [XI,YI] = MESHGRID(XI,YI)).
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD) where METHOD is one of
%       'linear'    - Triangle-based linear interpolation (default)
%       'cubic'     - Triangle-based cubic interpolation
%       'nearest'   - Nearest neighbor interpolation
%       'v4'        - MATLAB 4 griddata method
%   defines the type of surface fit to the data. The 'cubic' and 'v4'
%   methods produce smooth surfaces while 'linear' and 'nearest' have
%   discontinuities in the first and zero-th derivative respectively.  All
%   the methods except 'v4' are based on a Delaunay triangulation of the
%   data.
%   If METHOD is [], then the default 'linear' method will be used.
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD,OPTIONS) specifies a cell array 
%   of strings OPTIONS to be used as options in Qhull via DELAUNAYN. 
%   If OPTIONS is [], the default DELAUNAYN options will be used.
%   If OPTIONS is {''}, no options will be used, not even the default.
%
%   Example:
%      rand('seed',0)
%      x = rand(100,1)*4-2; y = rand(100,1)*4-2; z = x.*exp(-x.^2-y.^2);
%      ti = -2:.25:2; 
%      [xi,yi] = meshgrid(ti,ti);
%      zi = griddata(x,y,z,xi,yi);
%      mesh(xi,yi,zi), hold on, plot3(x,y,z,'o'), hold off
%
%   See also GRIDDATA3, GRIDDATAN, DELAUNAY, INTERP2, MESHGRID, DELAUNAYN.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.33.4.5 $  $Date: 2004/12/06 16:35:45 $

%function [zi] = func_fast_griddata(x,y,z,xi,yi,method,tri,td)
error(nargchk(8,8,nargin))

% Sort x and y so duplicate points can be averaged before passing to delaunay

%Need x,y and z to be column vectors
sz = numel(x);
x = reshape(x,sz,1);
y = reshape(y,sz,1);
z = reshape(z,sz,1);
sxyz = sortrows([x y z],[2 1]);
x = sxyz(:,1);
y = sxyz(:,2);
z = sxyz(:,3);

switch lower(method),
  case 'linear'
    zi = linear(x,y,z,xi,yi,tri,td);
  case 'cubic'
    zi = cubic(x,y,z,xi,yi,tri,td);
  case 'nearest'
    zi = nearest(x,y,z,xi,yi,tri,td);
  otherwise
    error('MATLAB:fast_griddata:UnknownMethod', 'Unknown method.');
end
  

%------------------------------------------------------------
function zi = linear(x,y,z,xi,yi,tri,td)
%LINEAR Triangle-based linear interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.


siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns


t = td ;
% Only keep the relevant triangles.
out = find(isnan(t));
if ~isempty(out), t(out) = ones(size(out)); end
tri = tri(t,:);

% Compute Barycentric coordinates (w).  P. 78 in Watson.
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
      (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
          (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
          (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
          (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
w(out,:) = zeros(length(out),3);

z = z(:).'; % Treat z as a row so that code below involving
            % z(tri) works even when tri is 1-by-3.
zi = sum(z(tri) .* w,2);

zi = reshape(zi,siz);

if ~isempty(out), zi(out) = NaN; end
%------------------------------------------------------------

%------------------------------------------------------------
function zi = cubic(x,y,z,xi,yi,tri,td)
%TRIANGLE Triangle-based cubic interpolation

%   Reference: T. Y. Yang, "Finite Element Structural Analysis",
%   Prentice Hall, 1986.  pp. 446-449.
%
%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

% Triangularize the data

% Find the nearest triangle (t)
t = td ;
zi = cubicmx(x,y,z,xi,yi,tri,t);
%------------------------------------------------------------

%------------------------------------------------------------
function zi = nearest(x,y,z,xi,yi,tri,td)
%NEAREST Triangle-based nearest neightbor interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these a columns
x = x(:); y = y(:); z = z(:); % Treat these as columns


% Find the nearest vertex
k = td ;
zi = k;
d = find(isfinite(k));
zi(d) = z(k(d));
zi = reshape(zi,siz);
%----------------------------------------------------------


