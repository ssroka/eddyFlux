clear;close all;clc
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
patch_str = 'Kur'; % 'GS'   'Kur'
% land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length
clrs = 'kkkkkkkkkk';
for yy = [3 7]
    
    switch yy
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
    
%     time = double([ncread(srcD,'time');ncread(srcJFM,'time')])*60*60; % hours
%     time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
    
%     lat = ncread(srcJFM,'latitude');
%     lon = ncread(srcJFM,'longitude');
    time.Year = yy+2000;
    switch patch_str
        case 'GS'
            % Gulf Stream
            lat_bnds = [25 45];
            lon_bnds = [275 305];
        case 'Kur'
            % Kurishio
            lat_bnds = [25 45];
            lon_bnds = [130 170];
    end
    
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(filename_data)
    
    g = var(SST_prime,0,3);
    V(yy) =nanmedian(g(:));
    figure(yy)
    contourf(lon(patch_lon),lat(patch_lat),nanmean(sshf_patch,3)',30,'k')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    
    
    avgSSHF(yy) = -nanmedian(sshf_patch(:));
    avgSLHF(yy) = -nanmedian(slhf_patch(:));
    T(yy) = time(end).Year;
    s_str = sprintf('SSHF [W/$$m^{-2}$$] (mean DJFM %d)',time(end).Year);
    L_str = sprintf('SLHF [W/$$m^{-2}$$] (mean DJFM %d)',time(end).Year);
    
    figure(100)
    subplot(3,1,1)
    plot(time(end).Year,nanmean(V(yy)),[clrs(yy) '*'],'markersize',10,'displayname',s_str)
    hold on
    title(sprintf('$$V$$(SST'') in time (median DJFM %d)',T(yy)),'interpreter','latex')
    subplot(3,1,2)
    plot(time(end).Year,avgSSHF(yy),[clrs(yy) '*'],'markersize',10,'displayname',s_str)
    hold on
    title(sprintf('ERA5 SSHF [W/$$m^{-2}$$] (median DJFM %d)',time(end).Year),'interpreter','latex')
    subplot(3,1,3)
    hold on
    plot(time(end).Year,avgSLHF(yy),[clrs(yy) '*'],'markersize',10,'displayname',L_str)
    title(sprintf('ERA5 SSHF [W/$$m^{-2}$$] (median DJFM %d)',time(end).Year),'interpreter','latex')

end
for i = 1:3
    subplot(3,1,i)
    xlabel('year')
% lh = legend('-dynamiclegend');
% set(lh,'interpreter','latex','location','best')
set(gca,'fontsize',20,'xlim',2000+[2,8])
set(gcf,'color','w','position',[ 24          35        1414         763])
end


%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/var_vs_flux','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/var_vs_flux_SST_prime','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/var_vs_flux_SST_prime_med','-dpng')


%}
    

    