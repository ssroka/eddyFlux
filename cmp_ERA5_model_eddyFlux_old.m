

subplot(3,2,1)
contourf(lon_box,lat_box,mean_sshf_ERA5_box)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')
clim_sshf = get(gca,'clim');

subplot(3,2,2)
contourf(lon_box,lat_box,mean_slhf_ERA5_box)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')
clim_slhf = get(gca,'clim');

subplot(3,2,3)
contourf(lon_box,lat_box,mean_model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15,'clim',clim_sshf)
title(sprintf('model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,4)
contourf(lon_box,lat_box,mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15,'clim',clim_slhf)
title(sprintf('model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,5)
contourf(lon_box,lat_box,mean_sshf_ERA5_box-mean_model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 - model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,6)
contourf(lon_box,lat_box,mean_slhf_ERA5_box - mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 -model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))

update_figure_paper_size()
if all(abs(abCD_factor-1)<1e-10)
    print(sprintf('%simgs/cmp_model_ERA5_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
else
    print(sprintf('%simgs/cmp_model_ERA5_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,strrep(num2str(abCD_factor),'.','_')),'-dpdf')
end
figure(1000)
subplot(2,2,1)
contourf(lon_box,lat_box,(mean_sshf_ERA5_box-mean_model_sshf)./mean_sshf_ERA5_box)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y) - model(x,y))/ERA5(x,y) SSHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(2,2,2)
contourf(lon_box,lat_box,(mean_slhf_ERA5_box - mean_model_slhf)./mean_slhf_ERA5_box)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y) -model(x,y))/ERA5(x,y) SLHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(2,2,3)
contourf(lon_box,lat_box,nanmean((sshf_model-sshf_ERA5_box)./sshf_ERA5_box,3)')
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y,t) - model(x,y,t))/ERA5(x,y,t) SSHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(2,2,4)
contourf(lon_box,lat_box,nanmean((slhf_model - slhf_ERA5_box)./slhf_ERA5_box,3)')
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y,t) -model(x,y,t))/ERA5(x,y,t) SLHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')

set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))

update_figure_paper_size()

if all(abs(abCD_factor-1)<1e-10)
    print(sprintf('%simgs/cmp_model_ERA5_rel_er_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
else
    print(sprintf('%simgs/cmp_model_ERA5_rel_er_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,strrep(num2str(abCD_factor),'.','_')),'-dpdf')
end



















