%% function varout = combtiles3D(fnm,ext,dirin,lenrec,IEEE,nx,ny,nz,halo,numl,IT,JT);
%  Maarten Buijsman, USM, 2024-1-15
%  Function extracts stacked 3D records and combines several tiles
%  In the case of 1 tile, buffers are also removed
% 
%  Input:
%  fnm: filename, e.g., fnm = [plat_190_blk_']
%  ext: extension, ext = '.BinF'
%  dirin: directory
%  lenrec: record length, lenrec = (nx+2*halo)*(ny+2*halo)*(nz)+2;  
%  IEEE: file format, IEEE = 'ieee-be'
%  nx,ny,nz,halo,numl: x, y, z tile dimensions, buffer, number of records
%  IT,JT: vectors of x and y tile numbers, as in e.g.,  
%  the loaded file is, fname = 
%  [fnm num2str(JT(jj)) '_' num2str(IT(ii)) ext] => plat_190_blk_24_37.BinF

%
%  Output:
%  varout: combined variable

function varout = combtiles3D(fnm,ext,dirin,lenrec,IEEE,nx,ny,nz,halo,numl,IT,JT);

% % test
% lenrec = lenrec2;
% halo   = 3;
% numl = MEIG;
% 
% IT = 37:38; JT = 22:24;
% % test

var3 = zeros(length(JT)*ny,length(IT)*nx,nz,numl);

js=1; je=ny;
for jj=1:length(JT)
    
    is  = 1; ie  = nx;
    inj = [js:je];
    
    for ii=1:length(IT)
        
        ini = [is:ie];
        
        fname = [fnm sprintf('%02d',JT(jj)) '_' sprintf('%02d',IT(ii)) ext]; % allows first digit to be a 0        
        
        %disp([dirin fname])
        
        fid    = fopen([dirin fname],'r',IEEE);
        
        if (fid==-1); disp(['fid = -1; read error in ' dirin fname]); end
                
        % extract all records
        for i=1:numl
            alldata = fread(fid,lenrec,'single');
%            length(alldata(2:end-1))/( (nx+halo*2)*(ny+halo*2) )
            var1 = permute(reshape(alldata(2:end-1),[nx+halo*2 ny+halo*2 nz]),[2 1 3]);
            varout(inj,ini,:,i) = var1(halo+1:ny+halo,halo+1:nx+halo,:);
        end
        
        fclose(fid);
        is=is+nx; ie=ie+nx;
        
        %figure; pcolor(var3(:,:,1)); shading flat; caxis([-0.005 0.005])        
        %figure; plot(var3(10,:,1))
        
    end
    
    js=js+ny; 
    je=je+ny;
    
end
