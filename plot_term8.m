
clear;clc;close all

year = 2003;

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')

L = 200000; % m

month_str = 'DJFM';

%%


load(sprintf('term68_%d_%d',L/1000,year),...
    'term6','term8','term6_const','term8_const','To_prime','DT_prime',...
    'DT_CTRL','Umag','Dq_CTRL','Dq_prime')

dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
load(dataFile,'time','lat','lon','patch_lat','patch_lon')

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);
    
% plot term 8

%         term6(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*as.*SST_prime.*U_mag_CTRL.*DT_patch_CTRL;
%         term8(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*U_mag_CTRL.*DT_diff_prime;

%        term6(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*aL.*SST_prime.*U_mag_CTRL.*q_diff_CTRL;
%         term8(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_CTRL.*q_diff_prime;
        

figure(1)
ax = subplot(2,5,1);
[~,h] = contourf(lon_er,lat_er,nanmean(term6(:,:,:,1),3)');
title(sprintf('term 6 $$\\rho_a c_p^{air} C_D^s \\alpha_s T_o'' U \\Delta T $$ $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax) 

ax = subplot(2,5,2);
[~,h] = contourf(lon_er,lat_er,nanmean(term6_const(:,:,:,1),3)');
title(sprintf('$$\\rho_a c_p^{air} C_D^s \\alpha_s $$ = %4.2f ',mode(mode(mode(term6_const(:,:,:,1))))),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,3);
[~,h] = contourf(lon_er,lat_er,nanmean(To_prime,3)');
title(sprintf('$$T_o''$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,4);
[~,h] = contourf(lon_er,lat_er,nanmean(Umag,3)');
title(sprintf('U'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,5);
[~,h] = contourf(lon_er,lat_er,nanmean(DT_CTRL,3)');
title(sprintf('$$\\Delta T $$'),'interpreter','latex')
format_fig(h,ax)

% latent 
ax = subplot(2,5,6);
[~,h] = contourf(lon_er,lat_er,nanmean(term6(:,:,:,2),3)');
title(sprintf('term 6 $$\\rho_a L_v C_D^L \\alpha_L T_o'' U \\Delta q $$ $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,7);
[~,h] = contourf(lon_er,lat_er,nanmean(term6_const(:,:,:,2),3)');
title(sprintf('$$\\rho_a L_v C_D^L \\alpha_L $$ = %4.2f ',mode(mode(mode(term6_const(:,:,:,2))))),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,8);
[~,h] = contourf(lon_er,lat_er,nanmean(To_prime,3)');
title(sprintf('$$T_o''$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,9);
[~,h] = contourf(lon_er,lat_er,nanmean(Umag,3)');
title(sprintf('U'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,10);
[~,h] = contourf(lon_er,lat_er,nanmean(Dq_CTRL,3)');
title(sprintf('$$\\Delta q$$'),'interpreter','latex')
format_fig(h,ax)


set(gcf,'color','w','position',[61 221 1359 581],'NumberTitle','off','Name',num2str(year))


%         term6(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*as.*SST_prime.*U_mag_CTRL.*DT_patch_CTRL;
%         term8(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*U_mag_CTRL.*DT_diff_prime;

%        term6(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*aL.*SST_prime.*U_mag_CTRL.*q_diff_CTRL;
%         term8(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_CTRL.*q_diff_prime;
        

figure(2)

% sensible 

ax = subplot(2,5,1);
[~,h] = contourf(lon_er,lat_er,nanmean(term8(:,:,:,1),3)');
title(sprintf('term 8 $$\\rho_a c_p^{air} C_D^s U \\Delta T'' $$ $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax) 

ax = subplot(2,5,2);
[~,h] = contourf(lon_er,lat_er,nanmean(term8_const(:,:,:,1),3)');
title(sprintf('$$\\rho_a c_p^{air} C_D^s $$ = %4.2f ',mode(mode(mode(term8_const(:,:,:,1))))),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,3);
[~,h] = contourf(lon_er,lat_er,nanmean(Umag,3)');
title(sprintf('$$U$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,4);
[~,h] = contourf(lon_er,lat_er,nanmean(DT_prime,3)');
title(sprintf('$$\\Delta T''$$'),'interpreter','latex')
format_fig(h,ax)

% latent

ax = subplot(2,5,6);
[~,h] = contourf(lon_er,lat_er,nanmean(term8(:,:,:,2),3)');
title(sprintf('term 8 $$\\rho_a L_v C_D^L  U \\Delta q'' $$ $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,7);
[~,h] = contourf(lon_er,lat_er,nanmean(term8_const(:,:,:,2),3)');
title(sprintf('$$\\rho_a L_v C_D^L \\alpha_L $$ = %4.2f ',mode(mode(mode(term8_const(:,:,:,2))))),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,8);
[~,h] = contourf(lon_er,lat_er,nanmean(Umag,3)');
title(sprintf('U'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,5,9);
[~,h] = contourf(lon_er,lat_er,nanmean(Dq_prime,3)');
title(sprintf('$$\\Delta q''$$'),'interpreter','latex')
format_fig(h,ax)

set(gcf,'color','w','position',[61 221 1359 581],'NumberTitle','off','Name',num2str(year))

figure(1)
update_figure_paper_size()
print(sprintf('imgs/deconstruct_term6_L_%d_%d',L/1000,year),'-dpdf')
figure(2)
update_figure_paper_size()
print(sprintf('imgs/deconstruct_term8_L_%d_%d',L/1000,year),'-dpdf')

function [] = format_fig(h,plt_num,max_val,min_val)

set(h,'edgecolor','none')
set(gca,'ydir','normal','fontsize',15)
colorbar
xlabel('deg')
ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end




