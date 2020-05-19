plt_error_box = false;

t_range = 1:484;

year = 2007;

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

add_ssh_flag = true;

%% begin

load(sprintf('model_n_ERA_data_%d',year))
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year));


if plt_error_box
    
    err_box_lat = [32 38];
    err_box_lon = [140 160];
    
    err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
    err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
    
    model_sshf = nanmean(model_full_sshf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    model_slhf = nanmean(model_full_slhf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    
    model_sshf_no_eddy = nanmean(model_no_eddy_sshf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    model_slhf_no_eddy = nanmean(model_no_eddy_slhf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    
    ERA5_sshf = nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    ERA5_slhf = nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    
    ERA5_sshf_CTRL = nanmean(era_no_eddy_sshf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    ERA5_slhf_CTRL = nanmean(era_no_eddy_slhf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    
    lat_er = lat(patch_lat);
    lat_er = lat_er(err_box_bnds_lat);
    lon_er = lon(patch_lon);
    lon_er = lon_er(err_box_bnds_lon);
    
else
    
    model_sshf = nanmean(model_full_sshf(:,:,t_range),3)';
    model_slhf = nanmean(model_full_slhf(:,:,t_range),3)';
    
    model_sshf_no_eddy = nanmean(model_no_eddy_sshf(:,:,t_range),3)';
    model_slhf_no_eddy = nanmean(model_no_eddy_slhf(:,:,t_range),3)';
    
    ERA5_sshf = nanmean(sshf_patch(:,:,t_range),3)';
    ERA5_slhf = nanmean(slhf_patch(:,:,t_range),3)';
    
    ERA5_sshf_CTRL = nanmean(era_no_eddy_sshf(:,:,t_range),3)';
    ERA5_slhf_CTRL = nanmean(era_no_eddy_slhf(:,:,t_range),3)';
    
    lat_er = lat(patch_lat);
    lon_er = lon(patch_lon);
    
end

%%
ax = subplot(2,2,1);
[h] = contourf(lon_er,lat_er,ERA5_sshf - ERA5_sshf_CTRL);
max_diff = max(ERA5_sshf(:) - ERA5_sshf_CTRL(:));
min_diff = min(ERA5_sshf(:) - ERA5_sshf_CTRL(:));
title(sprintf('ERA5 diff SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax,max_diff,min_diff)

ax = subplot(2,2,3);
[h] = plot(nanmean(ERA5_slhf - ERA5_slhf_CTRL,2),lat_er,'linewidth',2);
title(sprintf('ERA5 diff SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax,max_diff,min_diff)

ax = subplot(2,2,2);
[h] = contourf(lon_er,lat_er,model_sshf-model_sshf_no_eddy);
max_diff = max(model_sshf(:) - model_sshf_no_eddy(:));
min_diff = min(model_sshf(:) - model_sshf_no_eddy(:));
title(sprintf('model diff SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax,max_diff,min_diff)

ax = subplot(2,2,4);
[h] = plot(nanmean(model_slhf - model_slhf_no_eddy,2),lat_er,'linewidth',2);
title(sprintf('model diff SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax,max_diff,min_diff)

set(gcf,'color','w','position',[ 232  1  1188  801],'NumberTitle','off','Name',num2str(year))

if add_ssh_flag
    
    for i = [1 2]
        subplot(2,2,i)
        hold on
        plot_SSH_contour;
    end
    
end

update_figure_paper_size()

if add_ssh_flag
    print(sprintf('imgs/model_zonal_integral_ssh_%d',year),'-dpdf')
else
    print(sprintf('imgs/model_zonal_integral_%d',year),'-dpdf')
end

function [] = format_fig(h,plt_num,max_val,min_val)

set(gca,'ydir','normal','fontsize',15)
colorbar
ylabel('deg')
xlabel('zonal mean')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end