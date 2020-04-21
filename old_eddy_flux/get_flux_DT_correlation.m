clear;close all;clc
% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'box'; % 'box' or 'point'
err_box_lat = [32 38];
err_box_lon = [130 160];
error_units_flag = 'flux'; % 'flux' or 'power'

patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

alpha_eq_flag = 'constant_imposed'; % 'constant' (calculated from mean of flux (LS with constant valued model)
% 'linear'   (calculated from mean of flux (LS with model linear in T')
% 'constant_imposed' (pick the same constant for a_s and a_L

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length



for year = [3 7]
    
    switch year
        case 3 % 2003
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
        case 4 % 2004
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2004_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2003_Dec_31d_6h.nc';
        case 5 % 2005
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2005_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2004_Dec_31d_6h.nc';
        case 6 % 2006
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2006_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2005_Dec_31d_6h.nc';
        case 7 % 2007
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2007_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2006_Dec_31d_6h.nc';
    end
    
    time = double([ncread(srcD,'time');ncread(srcJFM,'time')])*60*60; % hours
    time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
    
    lat = ncread(srcJFM,'latitude');
    lon = ncread(srcJFM,'longitude');
    
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
    
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(filename_data)
    N = numel(DT_patch(:,:,1));
    if ~exist('pass_num','var')
        R2 = zeros(length(time),2);
        % use this branch to determine whether this is the first or
        % second time this loop is being run
        pass_num = 1;
        year1 = time(end).Year;
    else
        pass_num = 2;
        year2 = time(end).Year;
    end
    for tt = 1:length(time)
        %         sshf_as_constant(:,:,tt) = rho_a.*c_p_air.*CD_ref.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*(DT_patch(:,:,tt));
        %         slhf_aL_constant(:,:,tt) = rho_a.*Lv.*CD_ref.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*((qo_patch-qa_patch));
        [X,Y]  = create_grid(lon(patch_lon),lat(patch_lat));
        X_dist = abs(X(end,end))-(X(end,1));
        Y_dist = abs(Y(1,1)-Y(end,1));
        Nx     = floor(X_dist/L);
        Ny     = floor(Y_dist/L);
        sshf_patch_CTRL =  smooth2a(sshf_patch(:,:,tt),Ny, Nx);
        sshf_patch_CTRL(lsm_patch>0) = NaN;
        
        slhf_patch_CTRL =  smooth2a(slhf_patch(:,:,tt),Ny, Nx);
        slhf_patch_CTRL(lsm_patch>0) = NaN;
        
        
        SHF = sshf_patch(:,:,tt)-sshf_patch_CTRL+slhf_patch(:,:,tt)-slhf_patch_CTRL;
                DT = DT_patch(:,:,tt);
%         DT = SST_prime(:,:,tt);
        
        [B,~,~,~,STATS] = regress(SHF(:),[ones(N,1) DT(:)]);
        R2(tt,pass_num) = STATS(1);
    end
    
end

plot(time,R2,'linewidth',2)
xlabel('time (not year specific)')
ylabel('$$R^2 $$ for linear fit','interpreter','latex')
title('$$R^2$$ coefficient between $$Q_{Total}$$ and $$T''$$','interpreter','latex')
legend([num2str(year1) ' mean: ' num2str(mean(R2(:,1)))],...
    [num2str(year2) ' mean: ' num2str(mean(R2(:,2)))])
set(gca,'fontsize',24)
set(gcf,'color','w')



%{
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/R2_07_03_SST_Prime','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/R2_07_03_DT','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/R2_07_03_SST_Prime','-dpng')


%}