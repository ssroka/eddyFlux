
% close all
figure(2)
if ~(exist('err_box_bnds_lat','var'))
    err_box_bnds_lat = true(size(lat(patch_lat),1),1);
    err_box_bnds_lon = true(size(lon(patch_lon),1),1);
end

tt = 1:length(time);
lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

subplot(3,4,1)
contourf(lon_er,lat_er,nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('SSHF Wm$$^{-2}$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
sshf_clim = get(gca,'clim');

% subplot(3,4,2)
% contourf(lon_er,lat_er,nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
% set(gca,'ydir','normal','fontsize',20)
% title('SLHF Wm$$^{-2}$$ ERA5','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar
% slhf_clim = get(gca,'clim');


subplot(3,4,3)
contourf(lon_er,lat_er,nanmean(U_mag(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('U $$m/s$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,4)
contourf(lon_er,lat_er,nanmean(DT_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('DT K ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,5)
contourf(lon_er,lat_er,nanmean(sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('SSHF Wm$$^{-2}$$ model','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',sshf_clim)

% subplot(3,4,6)
% contourf(lon_er,lat_er,nanmean(slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
% set(gca,'ydir','normal','fontsize',20)
% title('SLHF Wm$$^{-2}$$ model','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar
% set(gca,'clim',slhf_clim)

subplot(3,4,7)
contourf(lon(err_box_bnds_lon),lat(err_box_bnds_lat),nanmean(SST_prime(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('SST'' K','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar


subplot(3,4,8)
contourf(lon_er,lat_er,nanmean(qo_patch(err_box_bnds_lon,err_box_bnds_lat,tt)-qa_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('$$q_o-q_a$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,9)
contourf(lon_er,lat_er,nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,tt)-sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',20)
title('SSHF ERA5 - model','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

% subplot(3,4,10)
% contourf(lon_er,lat_er,nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,tt)-slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
% set(gca,'ydir','normal','fontsize',20)
% title('SLHF ERA5 - model','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar

set(gcf,'color','w','position',[ 60          72        1300         733])

