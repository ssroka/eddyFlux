


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

lat_er_box_TF = lat>=er_box(1,1) & lat<=er_box(1,2);
lon_er_box_TF = lon>=er_box(2,1) & lon<=er_box(2,2);

% given a 2D array that is is of the box_opt size, use this to index out
% the error box elements
lat_box_2_er_TF = lat_er_box_TF(lat_box_TF);
lon_box_2_er_TF = lon_er_box_TF(lon_box_TF);
