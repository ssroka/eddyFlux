clear;close all;clc
% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'box'; % 'box' or 'point'
err_box_lat = [32 38];
err_box_lon = [130 160];
error_units_flag = 'flux'; % 'flux' or 'power'

patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';


% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length

clrs = hsv(5);
for year = 2003
    
    switch year
        case 2003 % 2003
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
        case 2004 % 2004
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2004_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2003_Dec_31d_6h.nc';
        case 2005 % 2005
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2005_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2004_Dec_31d_6h.nc';
        case 2006 % 2006
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2006_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2005_Dec_31d_6h.nc';
        case 2007 % 2007
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2007_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2006_Dec_31d_6h.nc';
    end
    
    time = double([ncread(srcD,'time');ncread(srcJFM,'time')])*60*60; % hours
    time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
    
    if ~exist('pass_num','var')
        time_master = time;
        pass_num = 1;
    end
    
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
    
    patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));
    
    m_LP = zeros(length(time),1);
    m_HP = zeros(length(time),1);
    m_NF = zeros(length(time),1);
    R2_LP = zeros(length(time),1);
    R2_HP = zeros(length(time),1);
    R2_NF = zeros(length(time),1);
    
    for tt = 1:length(time)
        Ta_patch_CTRL= smooth2a(t2m_patch(:,:,tt),Ny, Nx);
        Ta_patch_CTRL(lsm_patch>0) = NaN;
        Ta_prime = t2m_patch(:,:,tt) - Ta_patch_CTRL;
        
        SST_patch_CTRL = SST_patch(:,:,tt) - SST_prime(:,:,tt);
        
        inds_sm = find(patch_lat_sm(patch_lat));
        
        plot_set(:,:,1) = Ta_patch_CTRL;
        plot_set(:,:,2) = SST_patch_CTRL;
        plot_set(:,:,3) = Ta_prime;
        plot_set(:,:,4) = SST_prime(:,:,tt);
        
        Ta_LP_vec = plot_set(:,inds_sm,1);
        Ta_LP_vec = Ta_LP_vec(:);
        
        To_LP_vec = plot_set(:,inds_sm,2);
        To_LP_vec = To_LP_vec(:);
        
        Ta_HP_vec = plot_set(:,inds_sm,3);
        Ta_HP_vec = Ta_HP_vec(:);
        
        To_HP_vec = plot_set(:,inds_sm,4);
        To_HP_vec = To_HP_vec(:);
        
        Ta_NF_vec = t2m_patch(:,inds_sm,tt);
        To_NF_vec = SST_patch(:,inds_sm,tt);
        
        plot_set_vec(:,1) = Ta_LP_vec;
        plot_set_vec(:,2) = To_LP_vec;
        plot_set_vec(:,3) = Ta_HP_vec;
        plot_set_vec(:,4) = To_HP_vec;
        
        %         figure(1)
        %        plot_instance_HP_LP;
        
        
        %         figure(2)
        %         plot_instance_scatter;
        
        [B,~,~,~,STATS] = regress(To_LP_vec,[ones(length(Ta_LP_vec),1) Ta_LP_vec]);
        m_LP(tt) = B(2);
        R2_LP(tt) = STATS(1);
        
        [B,~,~,~,STATS] = regress(To_HP_vec,[ones(length(Ta_HP_vec),1) Ta_HP_vec]);
        m_HP(tt) = B(2);
        R2_HP(tt) = STATS(1);
        
        [B,~,~,~,STATS] = regress(To_NF_vec,[ones(length(Ta_NF_vec),1) Ta_NF_vec]);
        m_NF(tt) = B(2);
        R2_NF(tt) = STATS(1);
        
    end
    
    
    % plot_instance_R2;

figure(1)
plot(time-years(year-3),R2_LP,'color',clrs(year-2,:),'linestyle','-','linewidth',2,'displayname',['LPF' num2str(time(end).Year)])
hold on
plot(time-years(year-3),R2_HP,'color',clrs(year-2,:),'marker','*','linestyle','none','linewidth',2,'displayname',['HPF' num2str(time(end).Year)])
plot(time-years(year-3),R2_NF,'color',clrs(year-2,:),'marker','s','linestyle','none','linewidth',2,'displayname',['HPF' num2str(time(end).Year)])
title('$$R^2$$ statistic for LS regression ','interpreter','latex')
end
datetick('x','dd-mmm')
lh = legend('-dynamiclegend');
set(lh,'interpreter','latex');
xlabel('time')

set(gca,'ydir','normal','fontsize',20)
set(gcf,'color','w','position',[ 84          50        1156         7329])


















