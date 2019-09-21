clear;close all;clc
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

alpha_eq_flag = 'linear'; % 'constant' or 'linear'

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length

for year = 3
    
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
    
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(filename_data)
    
    [LAT,LON] = meshgrid(lat(patch_lat),lon(patch_lon));
    sig_X = 20;
    sig_Y = 3;
    X0 = 153;
    Y0 = 37;
    W = repmat(exp(-(LON-X0).^2/(2*sig_X.^2)-(LAT-Y0).^2/(2*sig_Y.^2)),1,1,length(time));
      W = SST_prime/max(SST_prime(:));
% W = ones(size(LAT,1),size(LAT,2),length(time));
    switch alpha_eq_flag
        case 'constant'
            as(isinf(as)) = NaN;
            aL(isinf(aL)) = NaN;
            as_beta = nanmedian(as(:));
            aL_beta = nanmedian(aL(:));
            as_est =@(SST_prime) as_beta;
            aL_est =@(SST_prime) aL_beta;
        case 'linear'
            as_NaN = as;
            as_NaN(isinf(as)) = NaN;
            as_NaN = as_NaN(:);
            ids_as = isnan(as_NaN);
            SST_NaN = SST_prime(:);
            ids_SST = isnan(SST_NaN);
            as_NaN(ids_as|ids_SST) = [];
            SST_NaN(ids_as|ids_SST) = [];
            Ws = W;
            Ws(ids_as|ids_SST) = [];
            X = [ones(numel(as_NaN),1) SST_NaN];
            as_beta = ((X'.*repmat(Ws,2,1))*X)\((X'.*repmat(Ws,2,1))*as_NaN);
            as_est =@(SST_prime)  0*(as_beta(1) + as_beta(2)*SST_prime);
            clear Ws
            
            aL_NaN = aL;
            aL_NaN(isinf(aL)) = NaN;
            aL_NaN = aL_NaN(:);
            ids_aL = isnan(aL_NaN);
            SST_NaN = SST_prime(:);
            ids_SST = isnan(SST_NaN);
            aL_NaN(ids_aL|ids_SST) = [];
            SST_NaN(ids_aL|ids_SST) = [];
            WL = W;
            WL(ids_aL|ids_SST) = [];
            X = [ones(numel(aL_NaN),1) SST_NaN];
            aL_beta = ((X'.*repmat(WL,2,1))*X)\((X'.*repmat(WL,2,1))*as_NaN);
            aL_est =@(SST_prime) 0*( aL_beta(1) + aL_beta(2)*SST_prime);
            clear WL
    end
    
    sshf_as_constant = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_aL_constant = zeros(sum(patch_lon),sum(patch_lat),length(time));
    sshf_er = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_er = zeros(sum(patch_lon),sum(patch_lat),length(time));
    
    for tt = 1:length(time)
        
        % ---- qo ----
        qo_patch = SAM_qsatWater(SST_patch(:,:,tt), P0_patch(:,:,tt)) ;
        qo_patch_CTRL = SAM_qsatWater(SST_patch_CTRL, P0_patch(:,:,tt)) ;
        
        % ---- qa ----
        e_sat = SAM_psatWater(t2m_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat;
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch = r./(1+r);
        
        sshf_as_constant(:,:,tt) = rho_a.*c_p_air.*CD_ref.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*(DT_patch(:,:,tt));
        slhf_aL_constant(:,:,tt) = rho_a.*Lv.*CD_ref.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*((qo_patch-qa_patch));
        
        sshf_er(:,:,tt) = sshf_as_constant(:,:,tt) + sshf_patch(:,:,tt);
        slhf_er(:,:,tt) = slhf_aL_constant(:,:,tt) + slhf_patch(:,:,tt);
        
    end
    
end


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
    contourf(lon(patch_lon),lat(patch_lat),-sshf_patch(:,:,n)')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(sprintf('ERA5 SSHF [W/$$m^{-2}$$]'),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)

    subplot(2,3,5)
    contourf(lon(patch_lon),lat(patch_lat),-slhf_patch(:,:,n)')
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
%}

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






%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_linear_noW','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_linear_Gaus','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_linear_SST','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/W_Gaus','-dpng')



%}