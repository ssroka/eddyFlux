clear;close all;clc



% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'point'; % 'box_sum' 'box_mean' or 'point'
err_box_lat = [32 38];
err_box_lon = [140 160];
% err_box_lon = [130 160];
L = 500000;%m

patch_str = 'Kur'; % 'GS'   'Kur'

switch patch_str
    case 'GS'
        % Gulf Stream
        lat_bnds = [25 45];
        lon_bnds = [275 305];
    case 'Kur'
        % Kurishio
        lat_bnds = [25 45];
        lon_bnds = [130 170];
        lat_sm_bnds = [30 40];
end


year = 2007;

filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'lat','lon','patch_lat','patch_lon','SST_patch','SST_prime')

[X,Y] = create_grid(lon(patch_lon),lat(patch_lat));
X_dist = abs(X(end,end))-(X(end,1));
Y_dist = abs(Y(1,1)-Y(end,1));
Nx     = floor(X_dist/L);
Ny     = floor(Y_dist/L);

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));

tt = 1:size(SST_prime,3);
lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

To = nanmean(SST_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)';
% Top = nanmean(SST_prime(err_box_bnds_lon,err_box_bnds_lat,tt),3)';


Top = smooth2a(To,Nx,Ny);
To_bar = To - Top;


subplot(2,2,1)
contourf(lon_er,lat_er,To)
set(gca,'ydir','normal','fontsize',15)
title('$$T_o$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(2,2,2)
contourf(lon_er,lat_er,Top)
set(gca,'ydir','normal','fontsize',15)
title('$$T_o''$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

Top_sm = zeros(size(Top));

% subplot(2,2,3)
% for i = 1:size(SST_prime,3)
% Top_sm(:,:,i) = smooth2a(Top(:,:,i),Ny,Nx);
% end

contourf(lon_er,lat_er,Top)
set(gca,'ydir','normal','fontsize',15)
title('$$\overline{T_o}$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

contourf(lon_er,lat_er,To-To_bar-Top)
set(gca,'ydir','normal','fontsize',15)
title('$$diff$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar









