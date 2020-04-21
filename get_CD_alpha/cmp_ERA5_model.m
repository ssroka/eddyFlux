

ERA5_sshf = nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
ERA5_slhf = nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';

model_sshf = nanmean(sshf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
model_slhf = nanmean(slhf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';


subplot(3,2,1)
contourf(lon_er,lat_er,ERA5_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,2)
contourf(lon_er,lat_er,ERA5_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,3)
contourf(lon_er,lat_er,model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,4)
contourf(lon_er,lat_er,model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,5)
contourf(lon_er,lat_er,ERA5_sshf-model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 - model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,6)
contourf(lon_er,lat_er,ERA5_slhf - model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 -model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))




















