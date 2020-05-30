function [mean_er] = optimize_alpha_CD(alpha_CD,year,err_box_bnds_lat,err_box_bnds_lon,L)

%{


%}

filename = sprintf('Qs_QL_optimization_data_L_%d_%d',L/1000,year);
load(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')

as = alpha_CD(1);
aL = alpha_CD(2);

CD_s = alpha_CD(3);
CD_L = alpha_CD(4);


%     as_multiplier = rho_a.*c_p_air.*U_mag.*(SST_patch-t2m_patch);
%     aL_multiplier = rho_a.*Lv.*U_mag.*(qo_patch-qa_patch);


sshf_model = CD_s.*(1+as.*SST_prime).*as_multiplier;
slhf_model = CD_L.*(1+aL.*SST_prime).*aL_multiplier;



er = (sshf_model(err_box_bnds_lon,err_box_bnds_lat,:) - sshf_patch(err_box_bnds_lon,err_box_bnds_lat,:)).^2+...
    + (slhf_model(err_box_bnds_lon,err_box_bnds_lat,:) - slhf_patch(err_box_bnds_lon,err_box_bnds_lat,:)).^2;

mean_er = sqrt(nanmean(nanmean(nanmean(er))));




end



%{
% load(filename_mult,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')


% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'point'; % 'box_sum' 'box_mean' or 'point'

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
%}
