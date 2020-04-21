

clear;close all;clc
% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'box'; % 'box' or 'point'
err_box_lat = [32 38];
err_box_lon = [130 160];
error_units_flag = 'flux'; % 'flux' or 'power'

patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';
lsm = ncread(land_src,'lsm');

alpha_eq_flag = 'constant_imposed'; % 'constant' (calculated from mean of flux (LS with constant valued model)
% 'linear'   (calculated from mean of flux (LS with model linear in T')
% 'constant_imposed' (pick the same constant for a_s and a_L

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length

aviso = '/Volumes/SydneySroka_Remy/eddyFlux_data/zos_AVISO_L4_199210-201012.nc';

for year = 1:2
    
    switch year
        case 1 % 2003
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
        case 2 % 2007
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2007_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2006_Dec_31d_6h.nc';
    end
    
    time = double([ncread(srcD,'time');ncread(srcJFM,'time')])*60*60; % hours
    time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
    
    lat = ncread(srcJFM,'latitude');
    lon = ncread(srcJFM,'longitude');
    
    switch patch_str
        case 'GS'
            % Gulf Stream
            lat_bnds = [25 45];
            lon_bnds = [275 305];
        case 'Kur'
            % Kurishio
            lat_bnds = [25 45];
            lon_bnds = [130 170];
    end
   
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(filename_data)
    figure(year)
[c, h] = contourf(lon(patch_lon),lat(patch_lat),nanmean(sshf_patch+slhf_patch,3)');
            set(h,'linecolor','none')
            hold on


ylabel('lat [deg]','interpreter','latex')
xlabel('lon [deg]','interpreter','latex')
set(gca,'ydir','normal','fontsize',20)
colorbar
set(gcf,'color','w')
            
end


SSH_m = ncread(aviso,'zos');
lat_m = ncread(aviso,'lat');
lon_m = ncread(aviso,'lon');
t = ncread(aviso,'time')*24*3600;
t = datetime(t,'ConvertFrom','epochtime','epoch','1950-01-01');

    patch_lat_m = (lat_m>lat_bnds(1))&(lat_m<lat_bnds(2));
    patch_lon_m = (lon_m>lon_bnds(1))&(lon_m<lon_bnds(2));
    

SSH_data = '/Volumes/SydneySroka_Remy/eddyFlux_data/CCAR_recon_sea_level_20000103_20090627_v1.nc';
time_ssh = ncread(SSH_data,'time')*24*3600;
time_ssh = datetime(time_ssh,'ConvertFrom','epochtime','epoch','1900-01-01');
SSH_a = ncread(SSH_data,'ssha'); % cm
lat_a = ncread(SSH_data,'lat');
lon_a = ncread(SSH_data,'lon');

patch_lat_a = (lat_a>lat_bnds(1))&(lat_a<lat_bnds(2));
patch_lon_a = (lon_a>lon_bnds(1))&(lon_a<lon_bnds(2));

count = 0;
for i = 1:length(time_ssh)
    Dec_02 = (time_ssh.Year(i)==2002) && (time_ssh.Month(i) == 12);
    JFM_03 = (time_ssh.Year(i)==2003) && (time_ssh.Month(i) <= 3);
    if Dec_02 || JFM_03
        figure(3)
        count = count +1;
        % find the value of t that is closest
        [mt, it] = min(abs(time_ssh(i)-t));
        [X, Y] = meshgrid(lon_m(patch_lon_m),lat_m(patch_lat_m));
        [Xq, Yq] = meshgrid(lon_a(patch_lon_a), lat_a(patch_lat_a));
        ssh_tot = interp2(X,Y,SSH_m(patch_lon_m,patch_lat_m,it)',Xq,Yq)+SSH_a(patch_lon_a,patch_lat_a,i)'/100;
        contour(lon_a(patch_lon_a),lat_a(patch_lat_a),ssh_tot,[1 1],'k')
    end
    Dec_06 = (time_ssh.Year(i)==2006) && (time_ssh.Month(i) == 12);
    JFM_07 = (time_ssh.Year(i)==2007) && (time_ssh.Month(i) <= 3);
    if Dec_06 || JFM_07
        figure(7)
        [mt, it] = min(abs(time_ssh(i)-t));
        [X, Y] = meshgrid(lon_m(patch_lon_m),lat_m(patch_lat_m));
        [Xq, Yq] = meshgrid(lon_a(patch_lon_a), lat_a(patch_lat_a));
        ssh_tot = interp2(X,Y,SSH_m(patch_lon_m,patch_lat_m,it)',Xq,Yq)+SSH_a(patch_lon_a,patch_lat_a,i)'/100;
        contour(lon_a(patch_lon_a),lat_a(patch_lat_a),ssh_tot,[1 1],'k')
    end
end
count
t_str = sprintf('mean DJFM Sea Surface Flux %d, weekly SSH contours at 100cm',time(end).Year);
title(t_str,'interpreter','latex')

%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/SSH_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/SSH_07','-dpng')


%}


