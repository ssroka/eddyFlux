

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
    'lat','lon','patch_lat','patch_lon','time','slhf_patch','sshf_patch');

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

HF = sshf_patch+slhf_patch;

[~,h] = contourf(lon_er,lat_er,nanmean(HF,3)');
set(h,'edgecolor','none')
set(gca,'clim',[50 500])
