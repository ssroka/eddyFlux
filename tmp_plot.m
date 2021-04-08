lat_plot = lat(box_lat);
lon_plot = lon(box_lon);

[~,h] = contourf(lon_plot,lat_plot,...
    nanmean(h_diff,3)');

qL  = Lv.*(qo_patch(opt_patch_lon,opt_patch_lat,tt)-qa_patch(opt_patch_lon,opt_patch_lat,tt));
cT = c_p_air.*(DT_patch(opt_patch_lon,opt_patch_lat,tt));
subplot(1,2,1)
[~,h] = contourf(lon_plot,lat_plot,...
    nanmean(qL./1000,3)');
colorbar

subplot(1,2,2)
[~,h] = contourf(lon_plot,lat_plot,...
    nanmean(cT./1000,3)');
colorbar





