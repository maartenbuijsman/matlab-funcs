clear all

dir1='/home/mbui/data/projects/carmen/3_hlp/';



B21 = nc_TtxtConvt3([dir1,'6b211t_3.nc']);
B22 = nc_TtxtConvt3([dir1,'6b221t_3.nc']);
B23 = nc_TtxtConvt3([dir1,'6b231t_3.nc']);
B24 = nc_TtxtConvt3([dir1,'6b241t_3.nc']);
B25 = nc_TtxtConvt3([dir1,'6b251t_3.nc']);
B26 = nc_TtxtConvt3([dir1,'6b261t_3.nc']);
B27 = nc_TtxtConvt3([dir1,'6b271t_3.nc']);
B28 = nc_TtxtConvt3([dir1,'6b281t_3.nc']);
B29 = nc_TtxtConvt3([dir1,'6b291t_3.nc']);
B2 = [B21;B22;B23;B24;B25;B26;B27;B28;B29];
Yin = B2(:,3);  %depth
Xin = B2(:,4);  %time
Zin = B2(:,5);  %temp


%
figure
plot(Xin,Yin,'k.')

figure
plot(Zin,'k.')


num_grid = 247*24; %time steps close to matching common hourly sampling rate
[Xgrid,Ygrid,Zgrid] = norm_interp2_calc(Xin,Yin,Zin,num_grid);


eval(['save /home/mbui/data/projects/peacock/output/SCs_out/HOV_',str,'_T.mat']); disp(['HOV_',str,'_T.mat saved'])
load(['fname_',num2str(i)])