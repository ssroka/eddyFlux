er_s = (sshf_model_opt - ERA5_sshf)./ERA5_sshf;
er_L = (slhf_model_opt - ERA5_slhf)./ERA5_slhf;

ind_2_find = find(abs(er_L)>10);

[Iv,Jv,Kv] = ind2sub(size(ERA5_slhf),ind_2_find);

for i = 1:length(Iv)
    I = Iv(i);
    J = Jv(i);
    K = Kv(i);
    plot(sshf_model_opt(I,J,K),ERA5_sshf(I,J,K),'*')
    hold on
    plot(slhf_model_opt(I,J,K), ERA5_slhf(I,J,K),'o')
end


