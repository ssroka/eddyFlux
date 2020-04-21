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
% alpha_vec = [0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10];
CD_ref_vec = [1e-5 5e-5 1e-4 5e-4 1e-3 5e-3];
alpha_vec = [0 0.005 0.05 0.5 5 10];

% alpha_vec = [0 5 10 50 100];

% alpha_vec = [1 10];
% CD_ref_vec = [1e-5 5e-4 5e-3];

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

for ii_alpha = 1:length(alpha_vec)
    fprintf('alpha = %f\n',alpha_vec(ii_alpha))
    for ii_CD = 1:length(CD_ref_vec)
        fprintf('\tC_D = %f\n',CD_ref_vec(ii_CD))
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
            clear CD_ref
            
            CD_ref = CD_ref_vec(ii_CD);
            
            patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));
            
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
                    as_est =@(SST_prime) alpha_vec(ii_alpha);
                    aL_est =@(SST_prime) alpha_vec(ii_alpha);
                    
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
                    as_est =@(SST_prime)  (as_beta(1) + as_beta(2)*SST_prime);
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
                    aL_est =@(SST_prime) ( aL_beta(1) + aL_beta(2)*SST_prime);
                    clear WL
            end
            
            sshf_as_constant = zeros(sum(patch_lon),sum(patch_lat),length(time));
            slhf_aL_constant = zeros(sum(patch_lon),sum(patch_lat),length(time));
            sshf_er = zeros(sum(patch_lon),sum(patch_lat),length(time));
            slhf_er = zeros(sum(patch_lon),sum(patch_lat),length(time));
            
            if strcmp(error_units_flag,'power')
                [X,Y] = create_grid(lon(patch_lon),lat(patch_lat));
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
                
                sshf_as_constant(:,:,tt) = rho_a.*c_p_air.*CD_ref.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*(DT_patch(:,:,tt));
                slhf_aL_constant(:,:,tt) = rho_a.*Lv.*CD_ref.*(1+aL_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*((qo_patch-qa_patch));
                
                sshf_er(:,:,tt) = (sshf_as_constant(:,:,tt) - sshf_patch(:,:,tt)).*dx.*dy;
                slhf_er(:,:,tt) = (slhf_aL_constant(:,:,tt) - slhf_patch(:,:,tt)).*dx.*dy;
                
            end
            
            
            
            if as_est(0) == aL_est(1) % if the alpha's are the same constant
                if ~exist('pass_num','var')
                    % use this branch to determine whether this is the first or
                    % second time this loop is being run
                    pass_num = 1;
                    tot_flux_1 = sshf_as_constant + slhf_aL_constant;
                    year1 = time(end).Year;
                else
                    pass_num = 2;
                    tot_flux_2 = sshf_as_constant + slhf_aL_constant;
                    year2 = time(end).Year;
                    
                end
                figure(ii_alpha*10+ii_CD)
                subplot(3,2,2*pass_num-1)
                inds_sm = find(patch_lat_sm(patch_lat));
                contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(sshf_as_constant(:,inds_sm,:) + slhf_aL_constant(:,inds_sm,:),3)',30,'k')
                colorbar
                xlabel(' Deg lon ')
                ylabel(' Deg lat ')
                p_str = '%';
                t_str = sprintf('Total Flux [%s] (mean DJFM %d)\n$$\\alpha$$ = % 4.3f $$\\qquad C_D^*$$ = %2.2e',error_units,time(end).Year,as_est(0),CD_ref);
                title(t_str,'interpreter','latex')
                set(gca,'ydir','normal','fontsize',20)
                set(gcf,'color','w')
                
                subplot(3,2,2*pass_num)
                contourf(lon(patch_lon),lat(patch_lat),nanmean(sshf_as_constant + slhf_aL_constant,3)',30,'k')
                colorbar
                xlabel(' Deg lon ')
                ylabel(' Deg lat ')
                p_str = '%';
                t_str = sprintf('Total Flux [%s] (mean DJFM %d)\n$$\\alpha$$ = % 4.3f $$\\qquad C_D^*$$ = %2.2e',error_units,time(end).Year,as_est(0),CD_ref);
                title(t_str,'interpreter','latex')
                set(gca,'ydir','normal','fontsize',20)
                
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
                
                figure(ii_alpha*100+ii_CD)
                if pass_num == 1
                    
                    if strcmp(error_metric_flag,'point')
                        
                        er_1_sm = sqrt(nanmean((tot_flux_1(:,inds_sm,:)-sshf_patch(:,inds_sm,:)-slhf_patch(:,inds_sm,:).^2),3));
                        er_1_Lg = sqrt(nanmean(tot_flux_1(:,:,:)-sshf_patch(:,:,:)-slhf_patch(:,:,:).^2,3));
                        
                        mean_er_1_sm(ii_alpha,ii_CD) = nanmean(nanmean(er_1_sm));
                        mean_er_1_Lg(ii_alpha,ii_CD) = nanmean(nanmean(er_1_Lg));
                        
                        
                        patch_lat_sm = patch_lat_sm;
                        patch_lon_sm = patch_lon;
                        
                    elseif strcmp(error_metric_flag,'box')
                        err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
                        err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
                        
                        tot_flux_1_box = tot_flux_1(err_box_bnds_lon,err_box_bnds_lat,:);
                        sshf_patch_box = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
                        slhf_patch_box = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
                        
                        er_1_sm = (nanmean(nanmean(tot_flux_1_box))- nanmean(nanmean(sshf_patch_box)) - nanmean(nanmean(slhf_patch_box)));
                      er_1_sm_time_avg = nanmean(tot_flux_1_box,3) - nanmean(sshf_patch_box,3) - nanmean(slhf_patch_box,3);
                        er_1_Lg = sqrt(nanmean((tot_flux_1(:,:,:)-sshf_patch(:,:,:)-slhf_patch(:,:,:)).^2,3));
                        
                        mean_er_1_sm(ii_alpha,ii_CD) = sqrt(nanmean((er_1_sm).^2));
                        mean_er_1_Lg(ii_alpha,ii_CD) = nanmean(nanmean(er_1_Lg));
                        
                        patch_lat_sm = (lat>err_box_lat(1)) & (lat<err_box_lat(2));
                        patch_lon_sm = (lon>err_box_lon(1)) & (lon<err_box_lon(2));
                        
                    end
                    
                    subplot(2,2,2*pass_num-1)
                    contourf(lon(patch_lon_sm),(lat(patch_lat_sm)),er_1_sm_time_avg',30,'k')
                    colorbar
                    xlabel(' Deg lon ')
                    ylabel(' Deg lat ')
                    p_str = '%';
                    t_str = sprintf('$$Q_{Total}(\\alpha$$ = % 4.3f $$, C_D^* =$$ %2.2e$$)$$ -ERA5 \n[%s] (mean DJFM %d)',as_est(0),CD_ref,error_units,time(end).Year);
                    title(t_str,'interpreter','latex')
                    set(gca,'ydir','normal','fontsize',20)
                    set(gcf,'color','w')
                    
                    subplot(2,2,2*pass_num)
                    contourf(lon(patch_lon),lat(patch_lat),er_1_Lg',30,'k')
                    colorbar
                    xlabel(' Deg lon ')
                    ylabel(' Deg lat ')
                    p_str = '%';
                    t_str = sprintf('$$Q_{Total}(\\alpha$$ = % 4.3f $$, C_D^* =$$ %2.2e$$)$$ -ERA5 \n[%s] (mean DJFM %d)',as_est(0),CD_ref,error_units,time(end).Year);
                    title(t_str,'interpreter','latex')
                    set(gca,'ydir','normal','fontsize',20)
                    
                    
                else
                    if strcmp(error_metric_flag,'point')
                        
                        er_2_sm = sqrt(nanmean((tot_flux_2(:,inds_sm,:)-sshf_patch(:,inds_sm,:)-slhf_patch(:,inds_sm,:)).^2,3));
                        er_2_Lg = nanmean(tot_flux_2(:,:,:)-sshf_patch(:,:,:)-slhf_patch(:,:,:),3);
                        
                        mean_er_2_sm(ii_alpha,ii_CD) = nanmean(nanmean(er_2_sm));
                        mean_er_2_Lg(ii_alpha,ii_CD) = nanmean(nanmean(er_2_Lg));
                        
                        patch_lat_sm = patch_lat_sm;
                        patch_lon_sm = patch_lon;
                        
                    elseif strcmp(error_metric_flag,'box')
                        err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
                        err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
                        
                        tot_flux_2_box = tot_flux_2(err_box_bnds_lon,err_box_bnds_lat,:);
                        sshf_patch_box = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
                        slhf_patch_box = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
                        
                        er_2_sm_time_avg = nanmean(tot_flux_2_box,3) - nanmean(sshf_patch_box,3) - nanmean(slhf_patch_box,3);
                        er_2_sm = (nanmean(nanmean(tot_flux_2_box))- nanmean(nanmean(sshf_patch_box)) - nanmean(nanmean(slhf_patch_box)));

                        er_2_Lg = nanmean(tot_flux_2(:,:,:)-sshf_patch(:,:,:)-slhf_patch(:,:,:),3);
                        
                        mean_er_2_sm(ii_alpha,ii_CD) = sqrt(nanmean((er_2_sm).^2));
                        mean_er_2_Lg(ii_alpha,ii_CD) = nanmean(nanmean(er_2_Lg));
                        
                        patch_lat_sm = (lat>err_box_lat(1)) & (lat<err_box_lat(2));
                        patch_lon_sm = (lon>err_box_lon(1)) & (lon<err_box_lon(2));
                        
                    end
                    
                    subplot(2,2,2*pass_num-1)
                    contourf(lon(patch_lon_sm),lat(patch_lat_sm),er_2_sm_time_avg',30,'k')
                    colorbar
                    xlabel(' Deg lon ')
                    ylabel(' Deg lat ')
                    p_str = '%';
                    t_str = sprintf('$$Q_{Total}(\\alpha$$ = % 4.3f $$, C_D^* =$$ %2.2e$$)$$ -ERA5 \n[%s] (mean DJFM %d)',as_est(0),CD_ref,error_units,time(end).Year);
                    title(t_str,'interpreter','latex')
                    set(gca,'ydir','normal','fontsize',20)
                    set(gcf,'color','w')
                    
                    subplot(2,2,2*pass_num)
                    contourf(lon(patch_lon),lat(patch_lat),er_2_Lg',30,'k')
                    colorbar
                    xlabel(' Deg lon ')
                    ylabel(' Deg lat ')
                    p_str = '%';
                    t_str = sprintf('$$Q_{Total}(\\alpha$$ = % 4.3f $$, C_D^* =$$ %2.2e$$)$$ -ERA5 \n[%s] (mean DJFM %d)',as_est(0),CD_ref,error_units,time(end).Year);
                    title(t_str,'interpreter','latex')
                    set(gca,'ydir','normal','fontsize',20)
                    set(gcf,'color','w','position',[ 84          50        1156         7329])
                end
                close all
            end
            close all
        end
        close all
        clear pass_num
    end
end


if strcmp(error_metric_flag,'box')
    error_type = sprintf('box average (%2.1f %2.1f)x(%2.1f %2.1f)',...
        err_box_lat(1),err_box_lat(2), err_box_lon(1),err_box_lon(2));
else
    error_type = 'mean point-wise';
end



figure(1000)
subplot(1,2,1)
contourf(CD_ref_vec,alpha_vec,mean_er_1_sm)
set(gca,'xscale','log')
t_str = sprintf('%s\nerror between ERA5 and $$Q_{Total}$$ %s (%d)',error_type,error_units,year1);
title(t_str,'interpreter','latex')
ylabel('$$\alpha$$','interpreter','latex')
xlabel('$$C_D^*$$','interpreter','latex')
set(gca,'ydir','normal','fontsize',20)
colorbar
set(gcf,'color','w')
subplot(1,2,2)
contourf(CD_ref_vec,alpha_vec,mean_er_2_sm)
set(gca,'xscale','log')
t_str = sprintf('%s\nerror between ERA5 and $$Q_{Total}$$ %s (%d)',error_type,error_units,year2);
title(t_str,'interpreter','latex')
ylabel('$$\alpha$$','interpreter','latex')
xlabel('$$C_D^*$$','interpreter','latex')
set(gca,'ydir','normal','fontsize',20)
colorbar
set(gcf,'color','w','position',[114         376        1220         422])

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


