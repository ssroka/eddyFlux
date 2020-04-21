
% close all
%{


filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'DT_patch')


%}


figure(year)
if ~(exist('err_box_bnds_lat','var'))
    err_box_bnds_lat = true(size(lat(patch_lat),1),1);
    err_box_bnds_lon = true(size(lon(patch_lon),1),1);
end

if ~(exist('DT_patch','var'))
filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'DT_patch')
end

tt = 1:size(SST_prime,3);
lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

subplot(2,3,1)
contourf(lon_er,lat_er,nanmean(term1(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('$$\rho_ac_pC_D\overline{||u||}\overline{T_o-T_a}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
sshf_clim = get(gca,'clim');

subplot(2,3,2)
contourf(lon_er,lat_er,nanmean(term2(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('$$\rho_ac_pC_D\alpha_s\overline{||u||}\overline{T_o''T_o''}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
sshf_clim = get(gca,'clim');

subplot(2,3,3)
contourf(lon_er,lat_er,nanmean(term3(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('$$-\rho_ac_pC_D\alpha_s\overline{||u||}\overline{T_o''T_a''}$$','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
sshf_clim = get(gca,'clim');

subplot(2,3,5)
contourf(lon_er,lat_er,nanmean(term2(err_box_bnds_lon,err_box_bnds_lat,tt),3)'./nanmean(term1(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('($$\rho_ac_pC_D\alpha_s\overline{||u||}\overline{T_o''T_o''}$$)/($$\rho_ac_pC_D\overline{||u||}\overline{T_o-T_a}$$)','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'colorscale','log')
sshf_clim = get(gca,'clim');

subplot(2,3,6)
contourf(lon_er,lat_er,nanmean(term3(err_box_bnds_lon,err_box_bnds_lat,tt),3)'./nanmean(term2(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('($$-\rho_ac_pC_D\alpha_s\overline{||u||}\overline{T_o''T_a''}$$)/($$\rho_ac_pC_D\alpha_s\overline{||u||}\overline{T_o''T_o''}$$)','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'colorscale','log')
sshf_clim = get(gca,'clim');

set(gcf,'color','w','position',[ -46          43        1450         762])
set(gcf,'numbertitle','off','name',num2str(year))


img_loc = '/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/';
% print([img_loc 'model_' num2str(year) '_eddy_contrib'],'-dpng')
% savefig([img_loc 'model_' num2str(year) '_eddy_contrib'])
% 
savefig([img_loc 'model_' num2str(year) '_terms_my_filter'])


% savefig([img_loc 'model_' num2str(year) '_eddy_contrib_4param'])
% savefig([img_loc 'model_' num2str(year) '_eddy_contrib_4param_smUmag'])
% savefig([img_loc 'model_' num2str(year) '_eddy_contrib_4param_smUmag_smdiff'])


