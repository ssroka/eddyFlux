function [mean_er] = optimize_alpha_beta_CD(abCD,year,lat_patch_2_box_TF,lon_patch_2_box_TF,L,filter_type,box_num,model_str)
%{


%}

filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year);

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


er = (sshf_model - sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)).^2+...
    + (slhf_model - slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)).^2;

mean_er = sqrt(nanmean(nanmean(nanmean(er))));


end



