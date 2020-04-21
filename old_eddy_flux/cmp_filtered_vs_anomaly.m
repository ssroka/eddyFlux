
figure(year-2000)
if ~(exist('err_box_bnds_lat','var'))
    err_box_bnds_lat = true(size(lat(patch_lat),1),1);
    err_box_bnds_lon = true(size(lon(patch_lon),1),1);
end

if ~(exist('DT_patch','var'))
filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'DT_patch')
end

lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

T_bar = SST_patch(:,:,tt) - T_o_prime;
T_prime = T_o_prime;
T_bar_sm = smooth_mat(T_bar,M);
T_prime_sm = smooth_mat(T_prime,M);



subplot(2,2,1)
contourf(lon_er,lat_er,T_bar')
set(gca,'ydir','normal','fontsize',15)
title('$$\overline{T}$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(2,2,2)
contourf(lon_er,lat_er,T_prime')
set(gca,'ydir','normal','fontsize',15)
title('$$T_o''$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(2,2,3)
contourf(lon_er,lat_er,T_bar_sm')
set(gca,'ydir','normal','fontsize',15)
title('$$\overline{\overline{T}}$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(2,2,4)
contourf(lon_er,lat_er,T_prime_sm')
set(gca,'ydir','normal','fontsize',15)
title('$$\overline{T_o''}$$ $$[K]$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar









