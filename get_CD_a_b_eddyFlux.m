
options = optimoptions(@fmincon,'Display','iter','steptolerance',1e-5);

%% begin optimizing for CD and alpha

X = zeros(length(abCD0),1);
FFINAL = zeros(length(abCD0),1);

cc = zeros(2,length(year_vec));
rms_er = zeros(2,length(year_vec));

for i = 1:length(year_vec)
    
    if calc_CD_alpha_flag
        
        [abCD,FFINAL] = fmincon(@(abCd) optimize_alpha_beta_CD(abCd,year_vec(i),lat_patch_2_box_TF,lon_patch_2_box_TF,L,filter_type,box_num,model_str),abCD0,[],[],[],[],LB,UB,[],options);
        
        save(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%s_%d',con_str,filter_type,L/1000,box_num,model_str,year_vec(i)),'abCD','FFINAL','box_opt');
        fprintf('saved %d\n',year_vec(i))
        [abCD0 abCD]
        
    end
    if plot_model_ERA_err_flag
        
        year = year_vec(i);
        
        filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year);
                
        t_range = 1:size(sshf_patch,3);
        
        load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%s_%d',con_str,filter_type,L/1000,box_num,model_str,year_vec(i)),'abCD','FFINAL','box_opt');
        
        
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
        
        sshf_ERA5_box = sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,t_range);
        slhf_ERA5_box = slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,t_range);
        
        mean_sshf_ERA5_box = nanmean(sshf_ERA5_box,3)';
        mean_slhf_ERA5_box = nanmean(slhf_ERA5_box,3)';
        
        figure(year)
        cmp_ERA5_model_eddyFlux;
        %             [cc(1,i),cc(2,i),rms_er(1,i),rms_er(2,i),rel_er_sshf,rel_er_slhf] = cmp_ERA5_model_corr(...
        %                 sshf_model_opt,slhf_model_opt,ERA5_sshf,ERA5_slhf);
        
        
    end
    
    
end




% if plot_model_ERA_err_flag
%     figure(2)
%     plot(year_vec,cc(1,:),'linewidth',2,'displayname','$$c$$ sensible')
%     hold on
%     plot(year_vec,cc(2,:),'linewidth',2,'displayname','$$c$$ latent')
%     xlabel('year','interpreter','latex')
%     ylabel('correlation coefficient','interpreter','latex')
%     legend('interpreter','latex','location','best')
%     set(gca,'fontsize',15)
%     set(gcf,'color','w')
%     update_figure_paper_size()
%     print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_cc_%d_%s_box%d_%d',L/1000,filter_type,box_num,year),'-dpdf')
%     
%     figure(3)
%     plot(year_vec,rms_er(1,:),'linewidth',2,'displayname','rms error sensible')
%     hold on
%     plot(year_vec,rms_er(2,:),'linewidth',2,'displayname','rms error latent')
%     xlabel('year','interpreter','latex')
%     ylabel('rms error [$$W/m^2$$]','interpreter','latex')
%     legend('interpreter','latex','location','best')
%     set(gca,'fontsize',15)
%     set(gcf,'color','w')
%     update_figure_paper_size()
%     print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_rms_%d_%s_box%d_%d',L/1000,filter_type,box_num,year),'-dpdf')
%     
% end

