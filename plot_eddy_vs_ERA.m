plt_error_box = false;

t_range = 1:484;

year = 2007;

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

add_ssh_flag = true;

%% begin

load(sprintf('model_n_ERA_data_%d',year))
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year));

if plt_error_box
    
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


% model_sshf = nanmean(sshf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
% model_slhf = nanmean(slhf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';

ax(1) = subplot(3,4,1);
contourf(lon_er,lat_er,ERA5_sshf)
colormap(ax(1),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 full SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(2) = subplot(3,4,2);
contourf(lon_er,lat_er,ERA5_slhf)
colormap(ax(2),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 full SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(3) = subplot(3,4,3);
contourf(lon_er,lat_er,model_sshf)
colormap(ax(3),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model full SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(4) = subplot(3,4,4);
contourf(lon_er,lat_er,model_slhf)
colormap(ax(4),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model full SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(5) = subplot(3,4,5);
contourf(lon_er,lat_er,ERA5_sshf_CTRL)
colormap(ax(5),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 no eddy SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(6) = subplot(3,4,6);
contourf(lon_er,lat_er,ERA5_slhf_CTRL)
colormap(ax(6),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 no eddy SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(7) = subplot(3,4,7);
contourf(lon_er,lat_er,model_sshf_no_eddy)
colormap(ax(7),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('no eddy SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(8) = subplot(3,4,8);
contourf(lon_er,lat_er,model_slhf_no_eddy)
colormap(ax(8),parula)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('no eddy SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(9) = subplot(3,4,9);
contourf(lon_er,lat_er,ERA5_sshf - ERA5_sshf_CTRL)
colormap(ax(9),rwb_map([max(ERA5_sshf(:) - ERA5_sshf_CTRL(:)) 0 min(ERA5_sshf(:) - ERA5_sshf_CTRL(:))],100))
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 diff SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(10) = subplot(3,4,10);
contourf(lon_er,lat_er,ERA5_slhf - ERA5_slhf_CTRL)
colormap(ax(10),rwb_map([max(ERA5_slhf(:) - ERA5_slhf_CTRL(:)) 0 min(ERA5_slhf(:) - ERA5_slhf_CTRL(:))],100))
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 diff SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(11) = subplot(3,4,11);
contourf(lon_er,lat_er,model_sshf-model_sshf_no_eddy)
colormap(ax(11),rwb_map([max(model_sshf(:) - model_sshf_no_eddy(:)) 0 min(model_sshf(:) - model_sshf_no_eddy(:))],100))
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model diff SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

ax(12) = subplot(3,4,12);
contourf(lon_er,lat_er,model_slhf - model_slhf_no_eddy)
colormap(ax(12),rwb_map([max(model_slhf(:) - model_slhf_no_eddy(:)) 0 min(model_slhf(:) - model_slhf_no_eddy(:))],100))
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model diff SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

if add_ssh_flag
    
for i = 1:12
    subplot(3,4,i)
    hold on
    plot_SSH_contour;
end

end

set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))




















