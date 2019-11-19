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
year_vec = [3 7];
CD_ref_vec = [0.00169 0.001555];
alpha_vec =  [1.27045 2.6630];


mean_er_1_sm  = zeros(length(alpha_vec),length(CD_ref_vec));
mean_er_1_Lg = zeros(length(alpha_vec),length(CD_ref_vec));
mean_er_2_sm  = zeros(length(alpha_vec),length(CD_ref_vec));
mean_er_2_Lg = zeros(length(alpha_vec),length(CD_ref_vec));
% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length


if strcmp(error_units_flag,'flux')
    error_units = sprintf('$$Wm^{-2}$$');
else
    error_units = sprintf('$$W$$');
end


for ii_year = 1:length(year_vec)
    
    switch year_vec(ii_year)
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
    clear CD_ref
    
    CD_ref = CD_ref_vec(ii_year);
    
    patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));
    inds_sm = find(patch_lat_sm(patch_lat));

    [LAT,LON] = meshgrid(lat(patch_lat),lon(patch_lon));
    [X,Y] = create_grid(lat(patch_lat),lon(patch_lon));
    sig_X = 20;
    sig_Y = 3;
    X0 = 153;
    Y0 = 37;
    W = repmat(exp(-(LON-X0).^2/(2*sig_X.^2)-(LAT-Y0).^2/(2*sig_Y.^2)),1,1,length(time));
    W = SST_prime/max(SST_prime(:));
    % W = ones(size(LAT,1),size(LAT,2),length(time));
    switch alpha_eq_flag
        case 'constant_imposed'
            as_est =@(SST_prime) alpha_vec(ii_year);
            aL_est =@(SST_prime) alpha_vec(ii_year);
        case 'constant'
            disp('not coded')
        case 'linear'
            disp('not coded')
    end
    
    sshf_eddy = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_eddy = zeros(sum(patch_lon),sum(patch_lat),length(time));
    sshf_sm = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_sm = zeros(sum(patch_lon),sum(patch_lat),length(time));
    sshf_diff = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_diff = zeros(sum(patch_lon),sum(patch_lat),length(time));
    SST_smooth =  zeros(sum(patch_lon),sum(patch_lat),length(time));
    SST_prime_smooth  =  zeros(sum(patch_lon),sum(patch_lat),length(time));
    
    [X,Y] = create_grid(lon(patch_lon),lat(patch_lat));
    X_dist = abs(X(end,end))-(X(end,1));
    Y_dist = abs(Y(1,1)-Y(end,1));
    Nx     = floor(X_dist/L);
    Ny     = floor(Y_dist/L);
    
    if strcmp(error_units_flag,'power')
        X_mid = 0.5*(X(:,1:end-1)+X(:,2:end));
        dx= [X_mid(:,1) X_mid(:,2:end)-X_mid(:,1:end-1) X(:,end)-X_mid(:,end)];
        Y_mid = 0.5*(Y(1:end-1,:)+Y(2:end,:));
        dy= [Y_mid(1,:); Y_mid(1:end-1,:)-Y_mid(2:end,:); Y(end,:)-Y_mid(end,:)];
    else
        dx = 1;
        dy = 1;
    end

                
    for tt = 1:length(time)
        
        % ---- qo ----
        qo_patch = SAM_qsatWater(SST_patch(:,:,tt), P0_patch(:,:,tt)) ;
        qo_patch_CTRL = SAM_qsatWater(SST_patch_CTRL, P0_patch(:,:,tt)) ;
        
        % ---- qa ----
        e_sat = SAM_psatWater(t2m_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat;
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch = r./(1+r);
        
        SST_patch_CTRL =  repmat(nanmean(SST_patch(:,:,tt)),sum(patch_lon),1); %
        SST_patch_CTRL(lsm_patch>0) = NaN;
        SST_smooth(:,:,tt) = SST_patch_CTRL;
        
        SST_smooth_CTRL = smooth2a(SST_smooth(:,:,tt),Ny, Nx);
        SST_smooth_CTRL(lsm_patch>0) = NaN;
        
        SST_prime_smooth(:,:,tt) = SST_smooth(:,:,tt) -  SST_smooth_CTRL;
        
        % ---- qo ----
        qo_patch_noeddy = SAM_qsatWater(SST_smooth(:,:,tt), P0_patch(:,:,tt)) ;
        qo_patch_CTRL_noeddy = SAM_qsatWater(SST_smooth_CTRL, P0_patch(:,:,tt)) ;
        
        
        DT_no_eddy = SST_smooth(:,:,tt)-t2m_patch(:,:,tt);
        
        sshf_eddy(:,:,tt) = rho_a.*c_p_air.*CD_ref.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*(DT_patch(:,:,tt));
        slhf_eddy(:,:,tt) = rho_a.*Lv.*CD_ref.*(1+aL_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*((qo_patch-qa_patch));
        
        sshf_sm(:,:,tt) = rho_a.*c_p_air.*CD_ref.*(1+as_est(SST_prime_smooth(:,:,tt)).*SST_prime_smooth(:,:,tt)).*U_mag(:,:,tt).*(DT_no_eddy);
        slhf_sm(:,:,tt) = rho_a.*Lv.*CD_ref.*(1+aL_est(SST_prime_smooth(:,:,tt)).*SST_prime_smooth(:,:,tt)).*U_mag(:,:,tt).*((qo_patch_noeddy-qa_patch));
                
        sshf_diff(:,:,tt) = (sshf_eddy(:,:,tt) - sshf_sm(:,:,tt));
        slhf_diff(:,:,tt) = (slhf_eddy(:,:,tt) - slhf_sm(:,:,tt));
        
%       plot_SST_primes;
    end
    
            
        tot_flux_1 = sshf_as_constant + slhf_aL_constant;
        
        
   

        
end


set(gca,'ydir','normal','fontsize',20)
colorbar

%{
0.001
5

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_linear_noW','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_linear_Gaus','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_linear_SST','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/W_Gaus','-dpng')


print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_0_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1e-3_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5e-3_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1e-2_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5e-2_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1e-1_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5e-1_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1_03_07','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5_03_07','-dpng')


print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_0_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1e-3_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5e-3_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1e-2_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5e-2_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1e-1_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5e-1_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_1_03_07_er','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/const_imp_a_5_03_07_er','-dpng')

    print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/mean_er_fxn_alpha','-dpng')

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/contour_a_CD_point_flux','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/contour_a_CD_box_flux','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/contour_a_CD_point_power','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/contour_a_CD_box_power','-dpng')


%}


% plot the linear version of the constant alpha and do the variance versus
% flux calculation tonight, then tomorrow morning you can evaluate what
% else you might want to do including weighted least squares, or you could
% do weighted least squares tonight.


%{
    n = 100;
    figure

    subplot(2,3,[1 4])
    contourf(lon(patch_lon),lat(patch_lat),SST_prime(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('SST'' [W/$$m^{-2}$$]\n(%d-%d-%d)',time(end).Year,time(end).Month,time(end).Day),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,2)
    contourf(lon(patch_lon),lat(patch_lat),sshf_patch(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('ERA5 SSHF [W/$$m^{-2}$$]'),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,5)
    contourf(lon(patch_lon),lat(patch_lat),slhf_patch(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('ERA5 SLHF [W/$$m^{-2}$$]'),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')

    subplot(2,3,3)
    contourf(lon(patch_lon),lat(patch_lat),sshf_as_constant(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SSHF $$\\alpha_s = %2.2f$$',as_constant),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,6)
    contourf(lon(patch_lon),lat(patch_lat),slhf_aL_constant(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SLHF $$\\alpha_L = %2.2f$$',aL_constant),'interpreter','latex')

    set(gca,'ydir','normal','fontsize',20)
     set(gcf,'color','w','position',[ 24          35        1414         763])

%%%%%%%%%%%%%%%%%%%%%%%%%

    n = 200;
    figure

    subplot(2,3,1)
    contourf(lon(patch_lon),lat(patch_lat),-sshf_patch(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('ERA5 SSHF [W/$$m^{-2}$$]'),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,4)
    contourf(lon(patch_lon),lat(patch_lat),-slhf_patch(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('ERA5 SLHF [W/$$m^{-2}$$]'),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    set(gcf,'color','w')

    subplot(2,3,2)
    contourf(lon(patch_lon),lat(patch_lat),sshf_as_constant(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SSHF $$\\alpha_s = %2.2f$$',as_constant),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,5)
    contourf(lon(patch_lon),lat(patch_lat),slhf_aL_constant(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SLHF $$\\alpha_L = %2.2f$$',aL_constant),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,3)
    contourf(lon(patch_lon),lat(patch_lat),sshf_as_constant(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SSHF difference '),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,6)
    contourf(lon(patch_lon),lat(patch_lat),slhf_aL_constant(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SLHF difference '),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)


    set(gca,'ydir','normal','fontsize',20)
     set(gcf,'color','w','position',[ 24          35        1414         763])

%%%%%%%%%%%%%%%%%%%%%%%%%


    figure

   switch alpha_eq_flag
        case 'constant'
            as_str = sprintf('median($$\\alpha_s$$) = %2.2f',as_beta(1))
            aL_str = sprintf('median($$\\alpha_L$$) = %2.2f',aL_beta(1))
        case 'linear'
            as_str = sprintf('%2.2f + %2.2f $$SST''$$',as_beta(1),as_beta(2))
            aL_str = sprintf('%2.2f + %2.2f $$SST''$$',aL_beta(1),aL_beta(2))
    end

    subplot(2,3,1)
    contourf(lon(patch_lon),lat(patch_lat),-nanmedian(sshf_patch,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    t_str = sprintf('ERA5 SSHF [W/$$m^{-2}$$]\n(median DJFM %d)',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,4)
    contourf(lon(patch_lon),lat(patch_lat),-nanmedian(slhf_patch,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    t_str = sprintf('ERA5 SLHF [W/$$m^{-2}$$]\n(median DJFM %d)',time(end).Year);
    title(t_str,'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,2)
    contourf(lon(patch_lon),lat(patch_lat),nanmedian(sshf_as_constant,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SSHF $$\\alpha_s = %2.2f$$\n(median DJFM %d)',as_constant,time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,5)
    contourf(lon(patch_lon),lat(patch_lat),nanmedian(slhf_aL_constant,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SLHF $$\\alpha_L = %2.2f$$\n(median DJFM %d)',aL_constant,time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,3)
    contourf(lon(patch_lon),lat(patch_lat),nanmedian(sshf_er,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SSHF difference\n(median DJFM %d)',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,6)
    contourf(lon(patch_lon),lat(patch_lat),nanmedian(slhf_er,3)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf(' SLHF difference\n(median DJFM %d)',time(end).Year),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)


    set(gca,'ydir','normal','fontsize',20)
     set(gcf,'color','w','position',[ 24          35        1414         763])

%%%%%%%%%%%%%%%%%%%%%%%%%


figure

switch alpha_eq_flag
    case 'constant'
        as_str = sprintf('median($$\\alpha_s$$) = %2.2f',as_beta(1));
        aL_str = sprintf('median($$\\alpha_L$$) = %2.2f',aL_beta(1));
    case 'linear'
        as_str = sprintf('%2.2e + %2.2e $$SST''$$',as_beta(1),as_beta(2));
        aL_str = sprintf('%2.2e + %2.2e $$SST''$$',aL_beta(1),aL_beta(2));
end

subplot(2,3,1)
contourf(lon(patch_lon),lat(patch_lat),-nanmedian(sshf_patch,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
t_str = sprintf('ERA5 SSHF [W/$$m^{-2}$$]\n(median DJFM %d)',time(end).Year);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,4)
contourf(lon(patch_lon),lat(patch_lat),-nanmedian(slhf_patch,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
t_str = sprintf('ERA5 SLHF [W/$$m^{-2}$$]\n(median DJFM %d)',time(end).Year);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,2)
contourf(lon(patch_lon),lat(patch_lat),nanmedian(sshf_as_constant,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
title(sprintf(' SSHF $$\\alpha_s$$ = %s \n(median DJFM %d)',as_str,time(end).Year),'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,5)
contourf(lon(patch_lon),lat(patch_lat),nanmedian(slhf_aL_constant,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
title(sprintf(' SLHF $$\\alpha_L $$ = %s\n(median DJFM %d)',aL_str,time(end).Year),'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,3)
contourf(lon(patch_lon),lat(patch_lat),nanmedian(sshf_er,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
title(sprintf(' SSHF difference\n(median DJFM %d)',time(end).Year),'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,6)
contourf(lon(patch_lon),lat(patch_lat),nanmedian(slhf_er,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
title(sprintf(' SLHF difference\n(median DJFM %d)',time(end).Year),'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)


set(gca,'ydir','normal','fontsize',20)
set(gcf,'color','w','position',[ 24          35        1414         763])



%}


