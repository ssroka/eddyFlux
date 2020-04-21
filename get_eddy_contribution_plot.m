

% ERA5_sshf = nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
% ERA5_slhf = nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
% 
% model_sshf = nanmean(sshf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
% model_slhf = nanmean(slhf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';


subplot(3,2,1)
contourf(lon_er,lat_er,model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model full SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,2)
contourf(lon_er,lat_er,model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model full SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,3)
contourf(lon_er,lat_er,model_sshf_no_eddy)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('no eddy SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,4)
contourf(lon_er,lat_er,model_slhf_no_eddy)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('no eddy SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,5)
contourf(lon_er,lat_er,model_sshf-model_sshf_no_eddy)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('full - no eddy SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,6)
contourf(lon_er,lat_er,model_slhf - model_slhf_no_eddy)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('full - no eddy SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))




















