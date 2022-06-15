%% function varout = combtiles2D(fnm,ext,dirin,lenrec,IEEE,nx,ny,halo,numl,IT,JT);
%  Maarten Buijsman, USM, 2022-6-15
%  Function extracts stacked 2D records and combines several tiles
%  In the case of 1 tile, buffers are also removed
% 
%  Input:
%  fnm: filename, fnm = ['modAN1_' RUNNM '_M2_JGFI_']
%  ext: extension, ext = '.BinF'
%  dirin: directory
%  lenrec: record length, lenrec = (nx+2*halo)*(ny+2*halo)+2;  
%  IEEE: file format, IEEE = 'ieee-be'
%  nx,ny,halo,numl: x and y tile dimensions, buffer, number of records
%  IT,JT: vectors of x and y tile numbers
%  Output:
%  varout: combined variable

% test
lenrec = lenrec2;
halo   = 3;
numl = MEIG;

IT = 37:38; JT = 22:24;

var3 = zeros(length(JT)*ny,length(IT)*nx,numl);

js=1; je=ny;
for jj=1:length(JT)
    
    is  = 1; ie  = nx;
    inj = [js:je];
    
    for ii=1:length(IT)
        
        ini = [is:ie];
        
        fname = [fnm num2str(JT(jj)) '_' num2str(IT(ii)) ext];
        fid    = fopen([dirin fname],'r',IEEE);
                
        % extract all records
        for i=1:numl
            alldata = fread(fid,lenrec,'single');
            var1 = permute(reshape(alldata(2:end-1),[nx+halo*2 ny+halo*2]),[2 1]);
            var3(inj,ini,i) = var1(halo+1:ny+halo,halo+1:nx+halo);
        end
        
        fclose(fid);
        is=is+nx; ie=ie+nx;
        
        %figure; pcolor(var3(:,:,1)); shading flat; caxis([-0.005 0.005])        
        %figure; plot(var3(10,:,1))
        
    end
    
    js=js+ny; je=je+ny;
    
end