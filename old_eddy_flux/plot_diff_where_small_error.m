clear;close all;clc

year = 2003;
thresh = 10; % dif W/m^2

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'point'; % 'box_sum' 'box_mean' or 'point'
err_box_lat = [32 38];
err_box_lon = [140 160];

patch_str = 'Kur'; % 'GS'   'Kur'

switch patch_str
    case 'GS'
        % Gulf Stream
        lat_bnds = [25 45];
        lon_bnds = [275 305];
    case 'Kur'
        % Kurishio
        lat_bnds = [25 45];
        lon_bnds = [130 170];
        lat_sm_bnds = [30 40];
end

%%

load(sprintf('%d_model_eddy_min_no_eddy_results.mat',year),'patch_lon','patch_lat',...
    'patch_str','lat','lon','lat_bnds','lon_bnds','U_mag','DT_patch','Lv',...
    'qo_patch','qa_patch','SST_prime','sshf_patch','slhf_eddy','sshf_eddy','slhf_sm','sshf_sm',...
    'slhf_patch','sshf_diff','slhf_diff','time','SST_patch','P0_patch','t2m_patch','RH_patch')

patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));


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

model_sL = slhf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT);
era5_sL  = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,TT);

era5_min_model_L = era5_sL-model_sL;

model_ss = sshf_eddy(err_box_bnds_lon,err_box_bnds_lat,TT);
era5_ss  = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,TT);

era5_min_model_s = era5_ss-model_ss;


subplot(3,4,3)
contourf(lon_er,lat_er,nanmean(era5_min_model_s,3)')
set(gca,'ydir','normal','fontsize',15)
title('ERA5 - model SSHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,4)
contourf(lon_er,lat_er,nanmean(era5_min_model_L,3)')
set(gca,'ydir','normal','fontsize',15)
title('ERA5 - model SLHF Wm$$^{-2}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

inds_s = abs(nanmean(era5_min_model_s,3)')<thresh;
inds_L = abs(nanmean(era5_min_model_L,3)')<thresh;

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

sensible_diff = nanmean(sshf_diff(err_box_bnds_lon,err_box_bnds_lat,TT),3)';
latent_diff = nanmean(slhf_diff(err_box_bnds_lon,err_box_bnds_lat,TT),3)';

sensible_diff(~inds_s) = NaN;
latent_diff(~inds_L) = NaN;


subplot(3,4,7)
contourf(lon_er,lat_er,sensible_diff)
set(gca,'ydir','normal','fontsize',15)
sm_dom_avg = round(nanmean(nanmean(nanmean(sensible_diff)')),2);
title(sprintf('SSHF full - no eddy Wm$$^{-2}$$, avg: %3.2f',sm_dom_avg),'interpreter','latex')

xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,8)
contourf(lon_er,lat_er,latent_diff)
set(gca,'ydir','normal','fontsize',15)
sm_dom_avg = round(nanmean(nanmean(nanmean(latent_diff)')),2);
title(sprintf('SLHF full - no eddy $$Wm^{-2}$$, avg: %3.2f',sm_dom_avg),'interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,9)
contourf(lon_er,lat_er,nanmean(sshf_sm(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
sm_dom_avg = round(nanmean(nanmean(nanmean(sensible_diff)')),2);
title(sprintf('SSHF no eddy Wm$$^{-2}$$, avg: %3.2f',sm_dom_avg),'interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(3,4,10)
contourf(lon_er,lat_er,nanmean(slhf_sm(err_box_bnds_lon,err_box_bnds_lat,TT),3)')
set(gca,'ydir','normal','fontsize',15)
sm_dom_avg = round(nanmean(nanmean(nanmean(sensible_diff)')),2);
title(sprintf('SLHF no eddy Wm$$^{-2}$$, avg: %3.2f',sm_dom_avg),'interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar








set(gcf,'color','w','position',[ -46          43        1450         762])
set(gcf,'numbertitle','off','name',num2str(year))








