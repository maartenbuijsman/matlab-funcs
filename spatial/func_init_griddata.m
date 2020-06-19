function [tri,ts,ds] = func_init_griddata(x,y,xi,yi,options)
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

%function [tri,ts,ds] = func_init_griddata(x,y,xi,yi,options)
error(nargchk(4,5,nargin))

if nargin == 5
    if ~iscellstr(options)
        error('MATLAB:OptsNotStringCell',...
              'OPTIONS should be cell array of strings.');           
    end
    opt = options;
else
    opt = [];
end

% Sort x and y so duplicate points can be averaged before passing to delaunay

%Need x,y and z to be column vectors
sz = numel(x);
x = reshape(x,sz,1);
y = reshape(y,sz,1);
sxy = sortrows([x y],[2 1]);
x = sxy(:,1);
y = sxy(:,2);
myeps = max(max(x)-min(x),max(y)-min(y))*1/2*eps^(1/3);
ind = [0; ((abs(diff(y)) < myeps) & (abs(diff(x)) < myeps)); 0];

if sum(ind) > 0
  warning('MATLAB:func_init_griddata:DuplicateDataPoints',['Duplicate x-y data points ' ...
            'detected: using average of the z values.']);
end

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x  = x(:);  y  = y(:);  % Treat these as columns

% Triangularize the data
if isempty(opt)
    tri = delaunayn([x y]);
else
    tri = delaunayn([x y],opt);
end
    
if isempty(tri),
  warning('MATLAB:func_init_griddata:CannotTriangulate','Data cannot be triangulated.');
  return
end

% Find the nearest triangle (t)
ts = tsearch(x,y,tri,xi,yi);

% Find the nearest vertex
ds = dsearch(x,y,tri,xi,yi);

%----------------------------------------------------------

