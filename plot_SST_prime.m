
function [] = plot_SST_prime(L,filter_type,model_str,box_num,lon_box,lat_box,data_src,year)
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
    'lat','lon','patch_lat','patch_lon','time');

filename = ...
    sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',...
    L/1000,filter_type,box_num,model_str,year);
load(filename,'SST_prime');

% lat_er = lat(patch_lat);
% lon_er = lon(patch_lon);

[~,h] = contourf(lon_box,lat_box,nanmean(SST_prime,3)');
set(h,'edgecolor','none')
set(gca,'clim',[-1 1]*2)
end