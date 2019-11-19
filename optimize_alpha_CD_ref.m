function [mean_er_1_sm] = optimize_alpha_CD_ref(alpha_CD,year)

%
% 
% non-smooth [alpha;cd]
% alpha = 1.270451301440938
% CD_ref = 0.001688020331327

%}



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

alpha_eq_flag = 'constant_imposed'; % 'constant' (calculated from mean of flux (LS with constant valued model)
% 'linear'   (calculated from mean of flux (LS with model linear in T')
% 'constant_imposed' (pick the same constant for a_s and a_L

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length


if strcmp(error_units_flag,'flux')
    error_units = sprintf('$$Wm^{-2}$$');
else
    error_units = sprintf('$$W$$');
end


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


CD_ref = alpha_CD(2);

as_est =@(SST_prime) alpha_CD(1);
aL_est =@(SST_prime) alpha_CD(1);

patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));

[LAT,LON] = meshgrid(lat(patch_lat),lon(patch_lon));
[X,Y] = create_grid(lat(patch_lat),lon(patch_lon));
sig_X = 20;
sig_Y = 3;
X0 = 153;
Y0 = 37;

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
    
    tot_flux_1 = sshf_as_constant + slhf_aL_constant;

    
    if strcmp(error_metric_flag,'point')
        
        er_1_sm = sqrt(nanmean((tot_flux_1(:,inds_sm,:)-sshf_patch(:,inds_sm,:)-slhf_patch(:,inds_sm,:).^2),3));
        
        mean_er_1_sm = nanmean(nanmean(er_1_sm));
        
        
    elseif strcmp(error_metric_flag,'box')
        err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
        err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
        
        tot_flux_1_box = tot_flux_1(err_box_bnds_lon,err_box_bnds_lat,:);
        sshf_patch_box = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
        slhf_patch_box = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
        
        er_1_sm = (nanmean(nanmean(tot_flux_1_box))- nanmean(nanmean(sshf_patch_box)) - nanmean(nanmean(slhf_patch_box)));
        
        mean_er_1_sm = sqrt(nanmean((er_1_sm).^2));
    end
    
    
    
    
    
end
end


