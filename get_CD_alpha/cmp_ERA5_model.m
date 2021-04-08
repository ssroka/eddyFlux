

subplot(3,2,1)
contourf(lon_er,lat_er,mean_ERA5_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')
clim_sshf = get(gca,'clim');

subplot(3,2,2)
contourf(lon_er,lat_er,mean_ERA5_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')
clim_slhf = get(gca,'clim');



subplot(3,2,3)
contourf(lon_er,lat_er,mean_model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15,'clim',clim_sshf)
title(sprintf('model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,4)
contourf(lon_er,lat_er,mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15,'clim',clim_slhf)
title(sprintf('model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,5)
contourf(lon_er,lat_er,mean_ERA5_sshf-mean_model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 - model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,6)
contourf(lon_er,lat_er,mean_ERA5_slhf - mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 -model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))

update_figure_paper_size()
print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_%d_%s_%d',L/1000,filter_type,year),'-dpdf')


figure(1000)
subplot(2,2,1)
contourf(lon_er,lat_er,(mean_ERA5_sshf-mean_model_sshf)./mean_ERA5_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y) - model(x,y))/ERA5(x,y) SSHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(2,2,2)
contourf(lon_er,lat_er,(mean_ERA5_slhf - mean_model_slhf)./mean_ERA5_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y) -model(x,y))/ERA5(x,y) SLHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(2,2,3)
contourf(lon_er,lat_er,nanmean((sshf_model_opt-ERA5_sshf)./ERA5_sshf,3)')
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y,t) - model(x,y,t))/ERA5(x,y,t) SSHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(2,2,4)
contourf(lon_er,lat_er,nanmean((slhf_model_opt - ERA5_slhf)./ERA5_slhf,3)')
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y,t) -model(x,y,t))/ERA5(x,y,t) SLHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))

update_figure_paper_size()
print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_rel_er_%d_%s_%d',L/1000,filter_type,year),'-dpdf')




















