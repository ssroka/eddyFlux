
% close all
%{


filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'DT_patch')


%}


figure(2)
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

subplot(4,4,1)
contourf(lon_er,lat_er,nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SSHF Wm$$^{-2}$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
sshf_clim = get(gca,'clim');

subplot(4,4,2)
contourf(lon_er,lat_er,nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SLHF Wm$$^{-2}$$ ERA5','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
slhf_clim = get(gca,'clim');


% subplot(4,4,3)
% contourf(lon_er,lat_er,nanmean(U_mag(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
% set(gca,'ydir','normal','fontsize',15)
% title('U $$m/s$$ ERA5','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar
% 
% subplot(4,4,4)
% contourf(lon_er,lat_er,nanmean(DT_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
% set(gca,'ydir','normal','fontsize',15)
% title('DT K ERA5','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar

model_sL = slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,tt);
era5_sL  = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,tt);

model_ss = sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,tt);
era5_ss  = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,tt);

model_sL_timeAvg = nanmean(model_sL,3);
era5_sL_timeAvg  = nanmean(era5_sL,3);

model_ss_timeAvg = nanmean(model_ss,3);
era5_ss_timeAvg  = nanmean(era5_ss,3);

SST_prime_er_box = SST_prime(err_box_bnds_lon,err_box_bnds_lat,tt);
DT_er_box = DT_patch(err_box_bnds_lon,err_box_bnds_lat,tt);

SST_prime_er_box_timeAvg = nanmean(SST_prime_er_box,3);
DT_er_box_timeAvg = nanmean(DT_er_box,3);

[B,Bint,~,~,STATS] = regress(era5_sL(:),[ones(numel(SST_prime_er_box),1) SST_prime_er_box(:)]);
[B,Bint,~,~,STATS_timeAvg] = regress(era5_ss_timeAvg(:),[ones(numel(SST_prime_er_box_timeAvg),1) SST_prime_er_box_timeAvg(:)]);

subplot(4,4,3)
plot(SST_prime_er_box(:),era5_ss(:),'bo','displayname','all')
hold on
plot(SST_prime_er_box_timeAvg(:),era5_ss_timeAvg(:),'ro','displayname','time avg')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('$$R^2_{all}$$ = %4.3f  $$R^2_{avg}$$ = %4.3f',STATS(1),STATS_timeAvg(1)),'interpreter','latex')
xlabel('ERA5 SST''')
ylabel('ERA5 SSHF')
lh = legend('-dynamiclegend');
set(lh,'location','best')


[B,Bint,~,~,STATS] = regress(era5_sL(:),[ones(numel(SST_prime_er_box),1) SST_prime_er_box(:)]);
[B,Bint,~,~,STATS_timeAvg] = regress(era5_sL_timeAvg(:),[ones(numel(SST_prime_er_box_timeAvg),1) SST_prime_er_box_timeAvg(:)]);

subplot(4,4,4)
plot(SST_prime_er_box(:),era5_sL(:),'bo','displayname','all')
hold on
plot(SST_prime_er_box_timeAvg(:),era5_sL_timeAvg(:),'ro','displayname','time avg')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('$$R^2_{all}$$ = %4.3f  $$R^2_{avg}$$ = %4.3f',STATS(1),STATS_timeAvg(1)),'interpreter','latex')
xlabel('ERA5 SST''')
ylabel('ERA5 SLHF')
lh = legend('-dynamiclegend');
set(lh,'location','best')


subplot(4,4,5)
contourf(lon_er,lat_er,nanmean(sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SSHF Wm$$^{-2}$$ model','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',sshf_clim)

subplot(4,4,6)
contourf(lon_er,lat_er,nanmean(slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SLHF Wm$$^{-2}$$ model','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar
set(gca,'clim',slhf_clim)

[B,Bint,~,~,STATS] = regress(era5_ss(:),[ones(numel(DT_er_box),1) DT_er_box(:) DT_er_box(:).^2]);
[B,Bint,~,~,STATS_timeAvg] = regress(era5_ss_timeAvg(:),[ones(numel(DT_er_box_timeAvg),1) DT_er_box_timeAvg(:)]);

subplot(4,4,7)
plot(DT_er_box(:),era5_ss(:),'bo','displayname','all')
hold on
plot(DT_er_box_timeAvg(:),era5_ss_timeAvg(:),'ro','displayname','time avg')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('$$R^2_{all}$$ = %4.3f  $$R^2_{avg}$$ = %4.3f',STATS(1),STATS_timeAvg(1)),'interpreter','latex')
xlabel('ERA5 DT')
ylabel('ERA5 SSHF')
lh = legend('-dynamiclegend');
set(lh,'location','best')

[B,Bint,~,~,STATS] = regress(era5_sL(:),[ones(numel(DT_er_box),1) DT_er_box(:) DT_er_box(:).^2]);
[B,Bint,~,~,STATS_timeAvg] = regress(era5_sL_timeAvg(:),[ones(numel(DT_er_box_timeAvg),1) DT_er_box_timeAvg(:)]);

subplot(4,4,8)
plot(DT_er_box(:),era5_sL(:),'bo','displayname','all')
hold on
plot(DT_er_box_timeAvg(:),era5_sL_timeAvg(:),'ro','displayname','time avg')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('$$R^2_{all}$$ = %4.3f  $$R^2_{avg}$$ = %4.3f',STATS(1),STATS_timeAvg(1)),'interpreter','latex')

xlabel('ERA5 DT')
ylabel('ERA5 SLHF')
lh = legend('-dynamiclegend');
set(lh,'location','best')





% subplot(4,4,8)
% contourf(lon_er,lat_er,nanmean(qo_patch(err_box_bnds_lon,err_box_bnds_lat,tt)-qa_patch(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
% set(gca,'ydir','normal','fontsize',15)
% title('$$q_o-q_a$$ ERA5','interpreter','latex')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% colorbar



subplot(4,4,9)
contourf(lon_er,lat_er,nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,tt)-sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SSHF ERA5 - model','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(4,4,10)
contourf(lon_er,lat_er,nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,tt)-slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SLHF ERA5 - model','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

subplot(4,4,11)
contourf(lon(err_box_bnds_lon),lat(err_box_bnds_lat),nanmean(SST_prime(err_box_bnds_lon,err_box_bnds_lat,tt),3)')
set(gca,'ydir','normal','fontsize',15)
title('SST'' K','interpreter','latex')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
colorbar

[B,Bint,~,~,STATS] = regress(era5_sL(:),[ones(numel(model_sL),1) model_sL(:)]);
[B,Bint,~,~,STATS_timeAvg] = regress(era5_ss_timeAvg(:),[ones(numel(model_ss_timeAvg),1) model_ss_timeAvg(:)]);

subplot(4,4,13)
plot(model_ss(:),era5_ss(:),'bo','displayname','all')
hold on
plot(model_ss_timeAvg(:),era5_ss_timeAvg(:),'ro','displayname','time avg')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('SSHF\n $$R^2_{all}$$ = %4.3f  $$R^2_{avg}$$ = %4.3f',STATS(1),STATS_timeAvg(1)),'interpreter','latex')
xlabel('model')
ylabel('ERA5')
lh = legend('-dynamiclegend');
set(lh,'location','best')

[B,Bint,~,~,STATS] = regress(era5_ss(:),[ones(numel(model_ss),1) model_ss(:)]);
[B,Bint,~,~,STATS_timeAvg] = regress(era5_sL_timeAvg(:),[ones(numel(model_sL_timeAvg),1) model_sL_timeAvg(:)]);

subplot(4,4,14)
plot(model_sL(:),era5_sL(:),'bo','displayname','all')
hold on
plot(model_sL_timeAvg(:),era5_sL_timeAvg(:),'ro','displayname','time avg')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('SLHF\n $$R^2_{all}$$ = %4.3f  $$R^2_{avg}$$ = %4.3f',STATS(1),STATS_timeAvg(1)),'interpreter','latex')
xlabel('model')
ylabel('ERA5')
lh = legend('-dynamiclegend');
set(lh,'location','best')

subplot(4,4,15) % histogram
histogram(nanmean(era5_ss,3)','displayname','ERA 5')
hold on
histogram(nanmean(sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)','displayname','model')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('SSHF DJFM %d',year),'interpreter','latex')
xlabel('Wm$$^{-2}$$ per $$\approx$$ 750 km$$^2$$','interpreter','latex')
legend('-dynamiclegend')


subplot(4,4,16) % histogram
histogram(nanmean(era5_sL,3)','displayname','ERA 5')
hold on
histogram(nanmean(slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,tt),3)','displayname','model')
set(gca,'ydir','normal','fontsize',15)
title(sprintf('SLHF DJFM %d',year),'interpreter','latex')
xlabel('Wm$$^{-2}$$ per $$\approx$$ 750 km$$^2$$','interpreter','latex')
legend('-dynamiclegend')


set(gcf,'color','w','position',[ -46          43        1450         762])
set(gcf,'numbertitle','off','name',num2str(year))


img_loc = '/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/';
% print([img_loc 'cmp_era_model_' num2str(year) '_1CD'],'-dpng')
print([img_loc 'cmp_era_model_' num2str(year) '_4params'],'-dpng')
savefig([img_loc 'cmp_era_model_' num2str(year) '_4params'])



