clear;close all;clc
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = true;
filter_flag        = 'box'; % 'box' or 'zonal'
patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length
% wind-related
rho_a = 1.2;     % kg m^-3
c_p_air = 1000;  % J / kg / K
Lv = 2.26e6;     % J/kg

CD_ref = 1e-3;   % reference drag coefficient

lsm = ncread(land_src,'lsm');

for year = 3:7
    close all
    clearvars -except year CD_ref L Lv c_p_air filter_flag flag_load_new_vars land_src lsm patch_str rho_a
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
    end
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    lsm_patch = lsm(patch_lon,patch_lat);
    
    patch_mat = sprintf('%s_%d_lat_%d_%d_lon_%d_%d.mat',patch_str,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));

    if flag_load_new_vars
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [slhf_patch] = get_patch('slhf',srcD,srcJFM,1./(-1*60*60),patch_lon,patch_lat,lsm_patch,true);
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [sshf_patch] = get_patch('sshf',srcD,srcJFM,1./(-1*60*60),patch_lon,patch_lat,lsm_patch,true);
    % SST [K]
    [SST_patch] = get_patch('sst',srcD,srcJFM,1,patch_lon,patch_lat,lsm_patch,true);
    % t2m [K]
    [t2m_patch] = get_patch('t2m',srcD,srcJFM,1,patch_lon,patch_lat,lsm_patch,true);
    % d2m [K]
    [d2m_patch] = get_patch('d2m',srcD,srcJFM,1,patch_lon,patch_lat,lsm_patch,true);
    % sp [Pa]
    [P0_patch] = get_patch('sp',srcD,srcJFM,1,patch_lon,patch_lat,lsm_patch,true);
    % u10 [m/s]
    [u10_patch] = get_patch('u10',srcD,srcJFM,1,patch_lon,patch_lat,lsm_patch,true);
    % v10 [m/s]
    [v10_patch] = get_patch('v10',srcD,srcJFM,1,patch_lon,patch_lat,lsm_patch,true);
    
    save(patch_mat,'slhf_patch','sshf_patch','SST_patch','t2m_patch','d2m_patch',...
         'P0_patch','u10_patch','v10_patch','time','lsm_patch');
    else
        load(patch_mat)
    end    
        
    SST_prime = zeros(sum(patch_lon),sum(patch_lat),length(time));
    DT_patch  = zeros(sum(patch_lon),sum(patch_lat),length(time));
    U_mag     = zeros(sum(patch_lon),sum(patch_lat),length(time));
    RH_patch  = zeros(sum(patch_lon),sum(patch_lat),length(time));
    Q_s  = zeros(sum(patch_lon),sum(patch_lat),length(time));
    Q_L  = zeros(sum(patch_lon),sum(patch_lat),length(time));
    
    as = zeros(sum(patch_lon),sum(patch_lat),length(time));
    aL = zeros(sum(patch_lon),sum(patch_lat),length(time));
    
    
    for tt = 1:length(time)
        
        switch filter_flag
            case 'box'
                [X,Y]  = create_grid(lon(patch_lon),lat(patch_lat));
                X_dist = abs(X(end,end))-(X(end,1));
                Y_dist = abs(Y(1,1)-Y(end,1));
                Nx     = floor(X_dist/L);
                Ny     = floor(Y_dist/L);
                
                SST_patch_CTRL = smooth2a(SST_patch(:,:,tt),Ny, Nx);
                SST_patch_CTRL(lsm_patch>0) = NaN;
                
                sshf_patch_CTRL =  smooth2a(sshf_patch(:,:,tt),Ny, Nx);
                sshf_patch_CTRL(lsm_patch>0) = NaN;
                
                slhf_patch_CTRL =  smooth2a(slhf_patch(:,:,tt),Ny, Nx);
                slhf_patch_CTRL(lsm_patch>0) = NaN;
            case 'zonal'
                [X,Y]  = create_grid(lon(patch_lon),lat(patch_lat));
                X_dist = abs(X(end,end))-(X(end,1));
                Y_dist = abs(Y(1,1)-Y(end,1));
                Nx     = floor(X_dist/L)+mod(floor(X_dist/L),2)-1;
                Ny     = floor(Y_dist/L)+mod(floor(Y_dist/L),2)-1;
                
                SST_patch_CTRL =  repmat(nanmean(SST_patch(:,:,tt)),sum(patch_lon),1); %
                SST_patch_CTRL(lsm_patch>0) = NaN;
                
                sshf_patch_CTRL =  repmat(nanmean(sshf_patch(:,:,tt)),sum(patch_lon),1); %
                sshf_patch_CTRL(lsm_patch>0) = NaN;
                
                slhf_patch_CTRL =  repmat(nanmean(slhf_patch(:,:,tt)),sum(patch_lon),1); %
                s1hf_patch_CTRL(lsm_patch>0) = NaN;
        end
        
        SST_prime(:,:,tt) = SST_patch(:,:,tt) - SST_patch_CTRL;
        
        DT_patch(:,:,tt) = SST_patch(:,:,tt) - t2m_patch(:,:,tt);
        
        U_mag(:,:,tt) = sqrt(u10_patch(:,:,tt).^2+v10_patch(:,:,tt).^2);
        
        RH_patch(:,:,tt)=100.*(exp((17.625.*(d2m_patch(:,:,tt)+273.15))./(243.04+(d2m_patch(:,:,tt)+273.15)))./exp((17.625.*(t2m_patch(:,:,tt)+273.15))./(243.04+(t2m_patch(:,:,tt)+273.15))));
        
        % ---- qo ----
        qo_patch = SAM_qsatWater(SST_patch(:,:,tt), P0_patch(:,:,tt)) ;
        qo_patch_CTRL = SAM_qsatWater(SST_patch_CTRL, P0_patch(:,:,tt)) ;
        
        % ---- qa ----
        e_sat = SAM_psatWater(t2m_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat;
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch = r./(1+r);
        
        % ---- CD ----
        CD = CD_ref;
        Q_s(:,:,tt) = rho_a.*c_p_air.*CD.*U_mag(:,:,tt).*(DT_patch(:,:,tt));
        Q_L(:,:,tt) = rho_a.*Lv.*CD.*U_mag(:,:,tt).*(qo_patch-qa_patch);
        
        as(:,:,tt) = (1-sshf_patch(:,:,tt)./(rho_a.*c_p_air.*CD_ref.*U_mag(:,:,tt).*(DT_patch(:,:,tt))))./SST_prime(:,:,tt);
        aL(:,:,tt) = (1-slhf_patch(:,:,tt)./(rho_a.*Lv.*CD_ref.*U_mag(:,:,tt).*(qo_patch-qa_patch)))./SST_prime(:,:,tt);
        
    end
    
    smooth_a;
    CDs = CD_ref.*(1+nanmean(SST_prime,3).*as_sm);
    CDL = CD_ref.*(1+nanmean(SST_prime,3).*aL_sm);
    filename_CD = sprintf('CD_%d_lat_%d_%d_lon_%d_%d.mat',time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    save(filename_CD,'CDs','CDL')
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    save(filename_data)

end



   %{
    contourf(X(patch_lat,patch_lon)./1000,Y(patch_lat,patch_lon)./1000,slhf_patch',10)
    colorbar
    title(['SST [K] (DJFM)  ' datestr(time(1),'yyyy')])
    set(gca,'ydir','reverse','xdir','reverse','fontsize',24)
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(SST_patch,3)',30,'k')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SST [K] (mean DJFM %d)',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(SST_prime,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SST'' [K] (mean DJFM %d)',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
        
%%% 
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmedian(SST_prime,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    t_str = sprintf('SST'' [K] (median DJFM %d)\n with %s SST'' ',time(end).Year,desc_str);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
%%%

    figure
    imagesc(lon(patch_lon),lat(patch_lat),nanmean(DT_patch,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' $$T_o-T_a$$ [K] (mean DJFM %d)',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(RH_patch,3)',30,'k')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    p_str = '%';
    t_str = sprintf('Relative Humidity $$[\\%s]$$ (mean DJFM %d)',char(37),time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(U_mag,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' $$\\| \\vec{u} \\|$$ [m/s] (mean DJFM %d)',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(slhf_patch,3)',30,'k')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    p_str = '%';
    t_str = sprintf('Surface Latent Heat Flux [W/$$m^{-2}$$] \n (mean DJFM %d) directly from ERA5',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(Q_s,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' $$Q_s$$ [W/$$m^{-2}$$] \n(mean DJFM %d) calculated with CD = 1E-3',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')

    figure
    contourf(lon(patch_lon),lat(patch_lat),nanmean(Q_L,3)',30,'k')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    p_str = '%';
    t_str = sprintf(' $$Q_L$$ [W/$$m^{-2}$$] \n(mean DJFM %d) calculated with CD = 1E-3',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),-nanmean(sshf_patch,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' Surface Sensible Heat Flux [W/$$m^{-2}$$] \n(mean DJFM %d) directly from ERA5',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')

    figure
    contourf(lon(patch_lon),lat(patch_lat),CDs')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' $$C_D^s$$ \n(from mean DJFM %d SST'')',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    contourf(lon(patch_lon),lat(patch_lat),CDL')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' $$C_D^L$$ \n(from mean DJFM %d SST'')',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    alpha_contour_levels = [-5 -1 1 5];
    contourf(lon(patch_lon),lat(patch_lat),-nanmean(as,3)',alpha_contour_levels)
    nColors = length(alpha_contour_levels); 
    colormap(parula(nColors));
    % Change this   ^ value as how many color reuired max is 64
    cbh = colorbar('eastoutside');
    cbh.Ticks = alpha_contour_levels;
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    t_str = sprintf('$$\\alpha_s$$ (mean DJFM %d) ',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'clim',[-10 10],'ydir','normal','fontsize',20)
    set(gcf,'color','w')
    
    figure
    alpha_contour_levels = [-5 -1 1 5];
    contourf(lon(patch_lon),lat(patch_lat),-nanmean(aL,3)',alpha_contour_levels)
    nColors = length(alpha_contour_levels); 
    colormap(parula(nColors));
    % Change this   ^ value as how many color reuired max is 64
    cbh = colorbar('eastoutside');
    cbh.Ticks = alpha_contour_levels;
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    t_str = sprintf('$$\\alpha_L$$ (mean DJFM %d) ',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'clim',[-10 10],'ydir','normal','fontsize',20)
    set(gcf,'color','w')

    %}



%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST_prime_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_box_SST_prime_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_DT_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_Umag_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_RH_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_sshf_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_sshf_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_Qs_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_Qs_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_slhf_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_slhf_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_QL_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_QL_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_as_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_aL_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_CDs_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_CDL_03','-dpng')


print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST_prime_zonal_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST_prime_zonal_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST_prime_box_03','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/era5_SST_prime_box_07','-dpng')


%}




