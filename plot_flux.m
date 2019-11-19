clear;close all;clc
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

% ----------- Plotting --------------
plt_flag(1) = false; % plot ERA5 SLHF and SSHF in 2x1 subplots

plt_flag(2) = true; % plot the ERA5 SLHF and SSHF in 2x2 grid where
% the left are the "small" domain (see lat_sm_bnds) and
% and the right are the "large" domain (see lat_bnds)
% and the top are the first year and the bottom are a
% second year.
%



if sum(double(plt_flag))>1
    error('I have not tested plotting more than one type')
end

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length
% wind-related
rho_a = 1.2;     % kg m^-3
c_p_air = 1000;  % J / kg / K
Lv = 2.26e6;     % J/kg

CD_ref = 1e-3;   % reference drag coefficient

lsm = ncread(land_src,'lsm');

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
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    patch_mat = sprintf('%s_%d_lat_%d_%d_lon_%d_%d.mat',patch_str,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    
    load(patch_mat)
    
    if plt_flag(1)
        figure(year)
        subplot(1,2,1)
        contourf(lon(patch_lon),lat(patch_lat),nanmean(slhf_patch,3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Surface Latent Heat Flux [W/$$m^{-2}$$] \n (mean DJFM %d) directly from ERA5',time(end).Year);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        
        subplot(1,2,2)
        contourf(lon(patch_lon),lat(patch_lat),nanmean(sshf_patch,3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Surface Sensible Heat Flux [W/$$m^{-2}$$] \n (mean DJFM %d) directly from ERA5',time(end).Year);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w','position',[84  333        1169         449])
        
    end
    
    if plt_flag(2)
        if ~exist('pass_num','var')
            % use this branch to determine whether this is the first or
            % second time this loop is being run
            pass_num = 1;
            tot_flux_1 = sshf_patch + slhf_patch;
            year1 = time(end).Year;
        else
            pass_num = 2;
            tot_flux_2 = sshf_patch + slhf_patch;
            year2 = time(end).Year;

        end
        figure(1)
        subplot(3,2,2*pass_num-1)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(slhf_patch(:,inds_sm,:)+sshf_patch(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Total ERA5 Flux [W/$$m^{-2}$$] (mean DJFM %d)',time(end).Year);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        subplot(3,2,2*pass_num)
        contourf(lon(patch_lon),lat(patch_lat),nanmean(slhf_patch+sshf_patch,3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Total ERA5 Flux [W/$$m^{-2}$$] (mean DJFM %d)',time(end).Year);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w','position',[84  333        1169         449])
        if pass_num == 2
        subplot(3,2,2*pass_num+1)
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(tot_flux_1(:,inds_sm,:)-tot_flux_2(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Difference (mean DJFM %d - %d)',year1,year2);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w','position',[84  333        1169         449])  
        subplot(3,2,2*pass_num+2)
        contourf(lon(patch_lon),lat(patch_lat),nanmean(tot_flux_1-tot_flux_2,3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Difference (mean DJFM %d - %d)',year1,year2);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w','position',[ 84          50        1156         7329])  
        end
    end
    
end


%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_mean_sshf_slhf_DJFM_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_mean_sshf_slhf_DJFM_04','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_mean_sshf_slhf_DJFM_05','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_mean_sshf_slhf_DJFM_06','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_mean_sshf_slhf_DJFM_07','-dpng')

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_total_SmLg_03_minus_07','-dpng')



%}