ccc
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

year = 2007;


load(sprintf('NCEP_data/NCEP_patch_data_%d.mat',year))
load(sprintf('NCEP_data/NCEP_patch_data_%d.mat',year))
whos

SST_patch

load('env_const.mat'); % load rho_a and c_p_air

CD = 0.0014;

m = size(SST_patch,1); % needs to be recalculated because of leap year
n = size(SST_patch,2); % needs to be recalculated because of leap year
p = size(SST_patch,3); % needs to be recalculated because of leap year
salinity = 34*ones(m,n,p);% ppt for Lv calculation

lsm = mean(U_mag,3)>5;

U_bar = U_mag;
Delta_h = SW_LatentHeat(SST_patch,'K',salinity,'ppt').*(qo_patch-qa_patch)+...
    c_p_air.*(SST_patch-t2m_patch);

beta_model_flux = mean(rho_a.*CD.*U_bar.*Delta_h ,3);
%%
NCEP_flux = mean((sshf_patch+slhf_patch),3);

subplot(2,5,1)
imagesc((beta_model_flux.*lsm)')
title('beta model with NCEP data and ERA5 CD')
colorbar
subplot(2,5,2)
imagesc((NCEP_flux.*lsm)')
title('NCEP flux [W/m^2]')
colorbar
subplot(2,5,3)
imagesc(mean((U_bar.*lsm),3)')
title('U')
colorbar
subplot(2,5,4)
imagesc(mean(Delta_h.*lsm,3)')
title('\Delta h [W/m^2]')
colorbar
subplot(2,5,5)
imagesc(mean(SST_patch.*lsm,3)')
colorbar
set(gca,'clim',[290 305])
title('SST [K]')
colorbar

%%
load(sprintf('ERA5_data/ERA5_patch_data_%d.mat',year))
whos

load('env_const.mat'); % load rho_a and c_p_air

CD = 0.0014;

m = size(SST_patch,1); % needs to be recalculated because of leap year
n = size(SST_patch,2); % needs to be recalculated because of leap year
p = size(SST_patch,3); % needs to be recalculated because of leap year
salinity = 34*ones(m,n,p);% ppt for Lv calculation


U_bar = U_mag;
Delta_h = SW_LatentHeat(SST_patch,'K',salinity,'ppt').*(qo_patch-qa_patch)+...
    c_p_air.*(SST_patch-t2m_patch);

beta_model_flux = mean(rho_a.*CD.*U_bar.*Delta_h ,3);
%%
ERA5_flux = mean(sshf_patch+slhf_patch,3);

subplot(2,5,6)
imagesc(beta_model_flux')
title('beta model with ERA5')
colorbar
subplot(2,5,7)
imagesc(ERA5_flux')
title('ERA5 flux [W/m^2]')
colorbar
subplot(2,5,8)
imagesc(mean(U_bar,3)')
title('U')
colorbar
subplot(2,5,9)
imagesc(mean(Delta_h,3)')
title('\Delta h [W/m^2]')
colorbar
subplot(2,5,10)
imagesc(mean(SST_patch,3)')
colorbar
set(gca,'clim',[290 305])
title('SST [K]')
colorbar
set(gcf,'position',[7         377        1428         428])
update_figure_paper_size()
print(sprintf('NCEP_as_filtered_data_%d',year),'-dpdf')













