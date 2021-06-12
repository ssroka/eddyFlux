ccc

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

year = 2003;


load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
    'SST_patch','lat','lon','patch_lat','patch_lon',...
    'slhf_patch','sshf_patch');

subplot(2,2,1)
[~,h] = contourf(lon(patch_lon),lat(patch_lat),SST_patch(:,:,1)');

hold on

% draw box
box_opt = [30 44.5; 148 169];

y_vec = [box_opt(1,1) box_opt(1,1) box_opt(1,2) box_opt(1,2) box_opt(1,1)];
x_vec = [box_opt(2,1) box_opt(2,2) box_opt(2,2) box_opt(2,1) box_opt(2,1)];

plot(x_vec,y_vec,'r','linewidth',3)


box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

lat_er = lat(box_lat);
lon_er = lon(box_lon);

% to index out of *_patch fields
inds_lat = box_lat(patch_lat);
inds_lon = box_lon(patch_lon);


lat_plot = lat(box_lat);
lon_plot = lon(box_lon);
subplot(2,2,2)
[~,h] = contourf(lon_plot,lat_plot,SST_patch(inds_lon,inds_lat,1)');

nans = any(isnan(SST_patch(inds_lon,inds_lat,1)),'all');
title(num2str(nans))

subplot(2,2,3)
[~,h] = contourf(lon(patch_lon),lat(patch_lat),SST_patch(:,:,1)');
hold on
% draw box
box_opt = [30 41.5; 142.5 169];

y_vec = [box_opt(1,1) box_opt(1,1) box_opt(1,2) box_opt(1,2) box_opt(1,1)];
x_vec = [box_opt(2,1) box_opt(2,2) box_opt(2,2) box_opt(2,1) box_opt(2,1)];

plot(x_vec,y_vec,'b','linewidth',3)


box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

lat_er = lat(box_lat);
lon_er = lon(box_lon);

% to index out of *_patch fields
inds_lat = box_lat(patch_lat);
inds_lon = box_lon(patch_lon);


lat_plot = lat(box_lat);
lon_plot = lon(box_lon);
subplot(2,2,4)
[~,h] = contourf(lon_plot,lat_plot,SST_patch(inds_lon,inds_lat,1)');

nans = any(isnan(SST_patch(inds_lon,inds_lat,1)),'all');
title(num2str(nans))



