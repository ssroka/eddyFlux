clear;close all;clc
addpath('~/Documents/MATLAB/util/')

patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

% ---------- Parameters --------------
% wind-related
rho_a = 1.2;     % kg m^-3
c_p_air = 1000;  % J / kg / K
Lv = 2.26e6;     % J/kg

CD_ref = 1e-3;   % reference drag coefficient
lsm = ncread(land_src,'lsm');

for year = 1
    
    switch year
        case 1 % 2003
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
        case 2 % 2007
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
    end
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    patch_mat = sprintf('%s_%d_lat_%d_%d_lon_%d_%d.mat',patch_str,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(patch_mat)
    load(sprintf('CD_%d_lat_%d_%d_lon_%d_%d.mat',time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2)),'CDs','CDL')

    Q_s_va  = zeros(sum(patch_lon),sum(patch_lat),length(time));
    Q_L_va  = zeros(sum(patch_lon),sum(patch_lat),length(time));
    
    
    lsm_patch = lsm(patch_lon,patch_lat);
    
    for tt = 1:length(time)
        
        SST_patch_CTRL =  repmat(nanmean(SST_patch(:,:,tt)),sum(patch_lon),1); %
        SST_patch_CTRL(lsm_patch>0) = NaN;
        SST_patch(:,:,tt) = SST_patch_CTRL;
        
        t2m_patch_CTRL =  repmat(nanmean(t2m_patch(:,:,tt)),sum(patch_lon),1); %
        t2m_patch_CTRL(lsm_patch>0) = NaN;
        t2m_patch(:,:,tt) = t2m_patch_CTRL;
                
        d2m_patch_CTRL =  repmat(nanmean(d2m_patch(:,:,tt)),sum(patch_lon),1); %
        d2m_patch_CTRL(lsm_patch>0) = NaN;
        d2m_patch(:,:,tt) = d2m_patch_CTRL;
        
        P0_patch_CTRL =  repmat(nanmean(P0_patch(:,:,tt)),sum(patch_lon),1); %
        P0_patch_CTRL(lsm_patch>0) = NaN;
        P0_patch(:,:,tt) = P0_patch_CTRL;
        
        u10_patch_CTRL =  repmat(nanmean(u10_patch(:,:,tt)),sum(patch_lon),1); %
        u10_patch_CTRL(lsm_patch>0) = NaN;
        u10_patch(:,:,tt) = u10_patch_CTRL;
        
        v10_patch_CTRL =  repmat(nanmean(v10_patch(:,:,tt)),sum(patch_lon),1); %
        v10_patch_CTRL(lsm_patch>0) = NaN;
        v10_patch(:,:,tt) = v10_patch_CTRL;

%         SST_prime(:,:,tt) = SST_patch(:,:,tt) - SST_patch_CTRL;
        
        DT_patch(:,:,tt) = SST_patch(:,:,tt) - t2m_patch(:,:,tt);
        
        U_mag(:,:,tt) = sqrt(u10_patch(:,:,tt).^2+v10_patch(:,:,tt).^2);
        
        RH_patch(:,:,tt)=100.*(exp((17.625.*(d2m_patch(:,:,tt)+273.15))./(243.04+(d2m_patch(:,:,tt)+273.15)))./exp((17.625.*(t2m_patch(:,:,tt)+273.15))./(243.04+(t2m_patch(:,:,tt)+273.15))));
        
        % ---- qo ----
        qo_patch = SAM_qsatWater(SST_patch(:,:,tt), P0_patch(:,:,tt));
        qo_patch_CTRL = SAM_qsatWater(SST_patch_CTRL, P0_patch(:,:,tt));
        
        % ---- qa ----
        e_sat = SAM_psatWater(SST_patch(:,:,tt)-DT_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat;
        r = 0.622 * e ./ max(e, P0_patch(:,:,tt)-e);
        qa_patch = r./(1+r);
        
        % ---- CD ----
        Q_s_va(:,:,tt) = rho_a.*c_p_air.*CDs.*U_mag(:,:,tt).*(DT_patch(:,:,tt));
        Q_L_va(:,:,tt) = rho_a.*Lv.*CDL.*U_mag(:,:,tt).*(qo_patch-qa_patch);
        

    end
    
 

end



%{

    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(Q_s_va,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('Vanished Anomaly $$Q_s$$ [W/$$m^{-2}$$] \n(mean DJFM %d) calculated with $$C_D^s$$',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')

    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(Q_L_va,3)',30,'k')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    p_str = '%';
    t_str = sprintf('Vanished Anomaly $$Q_L$$ [W/$$m^{-2}$$] \n(mean DJFM %d) calculated with $$C_D^L$$',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')



%}




%{
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_Qs_va_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_QL_va_03','-dpng')


print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_Qs_va_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_QL_va_07','-dpng')

%}




