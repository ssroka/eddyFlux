function [mean_er_1_sm] = optimize_alpha_CD_ref(alpha_CD,year)

%{
%
% non-smooth [alpha;cd]
% alpha = 1.270451301440938
% CD_ref = 0.001688020331327

alpha_CD
as0 = constant sensible heat coefficient
asL = linear sensible heat coefficient
aL0 = constant latent heat coefficient
aLL = linear latent heat coefficient

%}


% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'point'; % 'box_sum' 'box_mean' or 'point'
err_box_lat = [32 38];
err_box_lon = [140 160];
% err_box_lon = [130 160];

patch_str = 'Kur'; % 'GS'   'Kur'

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

filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_data,'lat','lon','patch_lat','patch_lon')

if length(alpha_CD) == 6
    CD_ref_s = alpha_CD(5);
    CD_ref_L = alpha_CD(6);
    
    as_est =@(SST_prime) alpha_CD(1)+alpha_CD(2)*abs(SST_prime);
    aL_est =@(SST_prime) alpha_CD(3)+alpha_CD(4)*abs(SST_prime);
    
elseif length(alpha_CD) == 5
    
    CD_ref_s = alpha_CD(5);
    CD_ref_L = alpha_CD(5);
    
    as_est =@(SST_prime) alpha_CD(1)+alpha_CD(2)*abs(SST_prime);
    aL_est =@(SST_prime) alpha_CD(3)+alpha_CD(4)*abs(SST_prime);
    
elseif length(alpha_CD) == 4
    
    CD_ref_s = alpha_CD(3);
    CD_ref_L = alpha_CD(4);
    
    as_est =@(SST_prime) alpha_CD(1);
    aL_est =@(SST_prime) alpha_CD(2);
    
end



patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));

% as_multiplier = rho_a.*c_p_air.*U_mag.*DT_patch;
% aL_multiplier = rho_a.*Lv.*U_mag.*(qo_patch-qa_patch);

filename_mult = sprintf('mult_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year-2000,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
load(filename_mult,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')

sshf_as_constant = CD_ref_s.*(1+as_est(SST_prime).*SST_prime).*as_multiplier;
slhf_aL_constant = CD_ref_L.*(1+aL_est(SST_prime).*SST_prime).*aL_multiplier;

tot_flux_1 = sshf_as_constant + slhf_aL_constant;


if strcmp(error_metric_flag,'point')
%     inds_sm = find(patch_lat_sm(patch_lat));
    err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
    err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
    
%     er_1_sm = (tot_flux_1(err_box_bnds_lon,err_box_bnds_lat,:)-sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:)-slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:)).^2;
    
    er = (sshf_as_constant(err_box_bnds_lon,err_box_bnds_lat,:) - sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:)).^2+...
            + (slhf_aL_constant(err_box_bnds_lon,err_box_bnds_lat,:) - slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:)).^2;
    
    mean_er_1_sm = sqrt(nanmean(nanmean(nanmean(er))));
    
    
elseif strcmp(error_metric_flag,'box_mean')
    err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
    err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
    
    tot_flux_1_box = tot_flux_1(err_box_bnds_lon,err_box_bnds_lat,:);
    sshf_patch_box = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
    slhf_patch_box = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
    
    er_1_sm = (nanmean(nanmean(tot_flux_1_box))- nanmean(nanmean(sshf_patch_box)) - nanmean(nanmean(slhf_patch_box)));
    
    mean_er_1_sm = sqrt(nanmean((er_1_sm).^2));
elseif strcmp(error_metric_flag,'box_sum')
    err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
    err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
    
    tot_flux_1_box = tot_flux_1(err_box_bnds_lon,err_box_bnds_lat,:);
    sshf_patch_box = sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
    slhf_patch_box = slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:);
    
    tot_flux_1_box(isnan(tot_flux_1_box)) = 0;
    sshf_patch_box(isnan(sshf_patch_box)) = 0;
    slhf_patch_box(isnan(slhf_patch_box)) = 0;
    
    er_1_sm = (sum(sum(tot_flux_1_box))- sum(sum(sshf_patch_box)) - sum(sum(slhf_patch_box)));
    
    mean_er_1_sm = sqrt(nanmean((er_1_sm).^2));
end


end


