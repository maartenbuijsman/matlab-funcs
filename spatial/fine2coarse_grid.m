function [outVar,outVar2] = fine2coarse_grid(lon,lat,LON,ddeg,gdepth,inVar,inVar2,inlon,inlat,inDepth,inAAA);

%% [outVar,outVar2] = fine2coarse_grid(lon,lat,LON,ddeg,gdepth,inVar,inVar2,inlon,inlat,inDepth,inAAA);
%  Maarten Buijsman, USM, 2019-6-6
%  average fine grid to coarse grid

% resolution of ingrid
inRes = inlon(500,2)-inlon(500,1);
nump = ceil(ddeg/inRes)*2;   %number of buffer points used
[ny,nx] = size(inlat);

outVar  = ones(size(LON))*NaN;
outVar2 = outVar;

for i=1:length(lon)-1
    
    % output to screen every 10 step
    if rem(i,10)==0
        disp(['i=' num2str(i) ' of ' num2str(length(lon)-1)]);
    end
    
    % polygon of the lon coordinates of the grid
    lonp = [lon(i) lon(i+1) lon(i+1) lon(i)   lon(i)];
    
    % ilat is the index of plat
    % it is advanced at the end of the loop
    ilat = 1;
    
    % loop over lat for constant lon
    for j=1:length(lat)-1
        
        %if ~isnan(gflux(j,i))
        if gdepth(j,i)>0 %then ocean; this retains more of the data
            % find plon(ilat,:) nearest to center of cell
            % cell center is ( lon(i) + lon(i+1) )/2
            % Iselx is index of longitude
            [d,Iselx] = min(abs(inlon(ilat,:) - (lon(i) + lon(i+1))/2));
            
            % find latitude corners of cell polygon
            latp = [lat(j) lat(j)   lat(j+1) lat(j+1) lat(j)];
            
            % speed up search
            % find plat(:,Iselx) nearest to center of cell
            % (lat(j) + lat(j+1))/2 is cell center
            % thus  Isely and Iselx are the indices on plon and plat closest to the cell center
            [d,Isely] = min(abs( inlat(:,Iselx) - (lat(j) + lat(j+1))/2 ));
            
            % select all indices in larger box marked by Isely and Iselx  +/- nump points
            II = max([Iselx-nump 1]):min([Iselx+nump nx]);
            JJ = max([Isely-nump 1]):min([Isely+nump ny]);
            %II,JJ, size(AAA);
            % find indices inside polygon, which marks the cell boundaries
            % isel is used for averaging
            isel = inpolygon(inlon(JJ,II),inlat(JJ,II),lonp,latp);
            
            % exlclude land and NaNs
            IN = find(isel==1 & inDepth(JJ,II)>0 & ~isnan(inVar(JJ,II)));
            numpa(j,i) = length(IN);
            
            % do averaging
            if ~isempty(IN)
                varTOTa2  = inVar(JJ,II);
                varTOTa22 = inVar2(JJ,II);
                AAA2     = inAAA(JJ,II);
                % area averaging
                outVar(j,i)  = sum(varTOTa2(IN).*AAA2(IN))/sum(AAA2(IN));
                outVar2(j,i) = sum(varTOTa22(IN).*AAA2(IN))/sum(AAA2(IN));
                
            end
            
            % test
            plot(lonp,latp,'r-')
            hold on
            plot(inlon(JJ,II),inlat(JJ,II),'k.')
            hold off
            pause
            
            % advance ilat
            % is this the problem????
            ilat = min([JJ(end) ny]);
            
        end
        
        %                 % advance ilat
        %                 ilat = min([JJ(end) ny]);
    end
end

% ------------------------------------------
