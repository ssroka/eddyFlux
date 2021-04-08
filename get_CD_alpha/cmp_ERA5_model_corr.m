function [cc_s,cc_L,rms_er_s,rms_er_L,rel_er_sshf,rel_er_slhf] = cmp_ERA5_model_corr(sshf_model,slhf_model,ERA5_sshf,ERA5_slhf);


cc_s = corr(sshf_model(:),ERA5_sshf(:));
cc_L = corr(slhf_model(:),ERA5_slhf(:));


rms_er_s = sqrt(mean((sshf_model(:)-ERA5_sshf(:)).^2));
rms_er_L = sqrt(mean((slhf_model(:)-ERA5_slhf(:)).^2));

rel_er_sshf = nanmean((sshf_model-ERA5_sshf)./ERA5_sshf,3)';
rel_er_slhf = nanmean((slhf_model-ERA5_slhf)./ERA5_slhf,3)';



