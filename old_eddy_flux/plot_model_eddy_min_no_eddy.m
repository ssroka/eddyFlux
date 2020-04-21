
% close all
%{


filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'DT_patch')


%}


figure(year-2000)
if ~(exist('err_box_bnds_lat','var'))
    err_box_bnds_lat = true(size(lat(patch_lat),1),1);
    err_box_bnds_lon = true(size(lon(patch_lon),1),1);
end

if ~(exist('DT_patch','var'))
filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'DT_patch')
end

TT = 1:size(SST_prime,3);
lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

subplot(3,4,1)
contourf(lon_er,lat_er,nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('ERA5 SSHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
sshf_clim = get(gca,'clim');

subplot(3,4,2)
contourf(lon_er,lat_er,nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('ERA5 SLHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
slhf_clim = get(gca,'clim');


% subplot(3,4,3)
% contourf(lon_er,lat_er,nanmean(U_mag(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
% set(gca,'ydir','normal','fontsize',15)
% title('U $$m/s$$ ERA5','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar
% 
% subplot(3,4,4)
% contourf(lon_er,lat_er,nanmean(DT_patch(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
% set(gca,'ydir','normal','fontsize',15)
% title('DT K ERA5','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar

model_sL = slhf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT);
era5_sL  = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,TT);

era5_min_model_L = era5_sL-model_sL;

model_ss = sshf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT);
era5_ss  = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,TT);

era5_min_model_s = era5_ss-model_ss;

% 
% model_sL_timeAvg = nanmean(model_sL,3);
% era5_sL_timeAvg  = nanmean(era5_sL,3);
% 
% model_ss_timeAvg = nanmean(model_ss,3);
% era5_ss_timeAvg  = nanmean(era5_ss,3);
% 
% SST_prime_er_box = SST_prime(err_box_bnds_lon,err_box_bnds_lat,TT);
% DT_er_box = DT_patch(err_box_bnds_lon,err_box_bnds_lat,TT);
% 
% SST_prime_er_box_timeAvg = nanmean(SST_prime_er_box,3);
% DT_er_box_timeAvg = nanmean(DT_er_box,3);

subplot(3,4,3)
contourf(lon_er,lat_er,nanmean(sshf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('model full SSHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',sshf_clim)

subplot(3,4,4)
contourf(lon_er,lat_er,nanmean(slhf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('model full SLHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',slhf_clim)

subplot(3,4,5)
contourf(lon_er,lat_er,nanmean(sshf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('model full SSHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',sshf_clim)

subplot(3,4,6)
contourf(lon_er,lat_er,nanmean(slhf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('model full SLHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',slhf_clim)


subplot(3,4,7)
contourf(lon_er,lat_er,nanmean(sshf_sm(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('model SSHF no eddy Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',sshf_clim)

subplot(3,4,8)
contourf(lon_er,lat_er,nanmean(slhf_sm(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
title('model SLHF no eddy Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',slhf_clim)

subplot(3,4,9)
contourf(lon_er,lat_er,nanmean(era5_min_model_s,3)')
set(gca,'ydir','normal','fontsize',15)
title('ERA5 - model SSHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,10)
contourf(lon_er,lat_er,nanmean(era5_min_model_L,3)')
set(gca,'ydir','normal','fontsize',15)
title('ERA5 - model SLHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,11)
contourf(lon_er,lat_er,nanmean(sshf_diff(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
sm_dom_avg = round(nanmean(nanmean(nanmean(sshf_diff(err_box_bnds_lon,err_box_bnds_lat,TT),3)')),2);
title(sprintf('SSHF full - no eddy Wm$$^{-2}$$, avg: %3.2f',sm_dom_avg),'interpreter','latex')

xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,12)
contourf(lon_er,lat_er,nanmean(slhf_diff(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
sm_dom_avg = round(nanmean(nanmean(nanmean(slhf_diff(err_box_bnds_lon,err_box_bnds_lat,TT),3)')),2);
title(sprintf('SLHF full - no eddy $$Wm^{-2}$$, avg: %3.2f',sm_dom_avg),'interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar



set(gcf,'color','w','position',[ -46          43        1450         762])
set(gcf,'numbertitle','off','name',num2str(year))


img_loc = '/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/';
% print([img_loc 'model_' num2str(year) '_eddy_contrib'],'-dpng')
% savefig([img_loc 'model_' num2str(year) '_eddy_contrib'])
% 

savefig([img_loc 'model_' num2str(year) '_eddy_contrib_impr_sm_my_filter'])
% savefig([img_loc 'model_' num2str(year) '_eddy_contrib_4param_smUmag'])
% savefig([img_loc 'model_' num2str(year) '_eddy_contrib_4param_smUmag_smdiff'])


