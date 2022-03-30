mean_rms = zeros(length(year_vec),1);

for i = 1:length(year_vec)
year = year_vec(i);

filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%s_%d',L/1000,filter_type,box_num,model_str,reanalysis_src,year);
load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%d_%s_%s_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,reanalysis_src,year_vec(i)),'abCD','FFINAL','box_opt');

abCD = abCD.*abCD_factor;

switch model_str
    case 'alpha'
        load(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')
        a = abCD(1);
        CD =  abCD(2);
        
        sshf_model = CD.*(1+a.*SST_prime).*as_multiplier;
        slhf_model = CD.*(1+a.*SST_prime).*aL_multiplier;
    case 'beta'
        load(filename,'bs_multiplier','bL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch')
        
        b = abCD(1);
        CD = abCD(2);
        sshf_model = CD.*(U_bar+b.*SST_prime).*bs_multiplier;
        slhf_model = CD.*(U_bar+b.*SST_prime).*bL_multiplier;
    case 'alphabeta'
        load(filename,'abs_multiplier','abL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch')
        
        a = abCD(1);
        b = abCD(2);
        CD = abCD(3);
        sshf_model = CD.*(1+a.*SST_prime).*(U_bar+b.*SST_prime).*abs_multiplier;
        slhf_model = CD.*(1+a.*SST_prime).*(U_bar+b.*SST_prime).*abL_multiplier;
end
mean_model_sshf = nanmean(sshf_model,3)';
mean_model_slhf = nanmean(slhf_model,3)';

sshf_ERA5_box = sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:);
slhf_ERA5_box = slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:);

mean_sshf_ERA5_box = nanmean(sshf_ERA5_box,3)';
mean_slhf_ERA5_box = nanmean(slhf_ERA5_box,3)';

mean_rms = sqrt(mean((mean_model_sshf(:)+mean_model_slhf(:)-...
    (mean_sshf_ERA5_box(:)+mean_slhf_ERA5_box(:))).^2));
save(sprintf('mean_rms_err_%s_%s_%d%s',model_str,reanalysis_src,year,abCD_fac_str),'mean_rms')
end


