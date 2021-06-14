


load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'lat','lon','patch_lat','patch_lon');

lat_box_TF = lat>=box_opt(1,1) & lat<=box_opt(1,2);
lon_box_TF = lon>=box_opt(2,1) & lon<=box_opt(2,2);

% to index out of *_patch fields
lat_patch_2_box_TF = lat_box_TF(patch_lat);
lon_patch_2_box_TF = lon_box_TF(patch_lon);

lat_box = lat(lat_box_TF);
lon_box = lon(lon_box_TF);

n = length(lat_box);
m = length(lon_box);
