

%% parameters for optimization




if alpha_pos_flag
    LB = [0;0; 0];
    UB = [Inf; Inf; Inf];
    con_str = 'cons_';
else
    LB = [-Inf; -Inf; 0];
    UB = [Inf; Inf; Inf];
    con_str = '';
end

options = optimoptions(@fmincon,'Display','iter','steptolerance',1e-5);

%% begin optimizing for CD and alpha

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'lat','lon','patch_lat','patch_lon');

file_for_box = load(sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,2003),'box_opt');

box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

% to index out of *_patch fields
opt_patch_lat = box_lat(patch_lat);
opt_patch_lon = box_lon(patch_lon);

prime_lat =  lat>=file_for_box.box_opt(1,1) & lat<=file_for_box.box_opt(1,2);
prime_lon = lon>=file_for_box.box_opt(2,1) & lon<=file_for_box.box_opt(2,2);

% to index out of *_prime fields
opt_prime_lat = box_lat(prime_lat);
opt_prime_lon = box_lon(prime_lon);

X = cell(length(year_vec),1);
FFINAL = cell(length(year_vec),1);

cc = zeros(2,length(year_vec));
rms_er = zeros(2,length(year_vec));

for i = 1:length(year_vec)
    
    if calc_CD_alpha_flag
        
        [x,ffinal] = fmincon(@(aCd) optimize_alpha_CD(aCd,year_vec(i),opt_prime_lat,opt_prime_lon,opt_patch_lat,opt_patch_lon,L,filter_type,box_num),aCd0,[],[],[],[],LB,UB,[],options);
        
        X{i} = x;
        FFINAL{i}  =ffinal;
        save(sprintf('opt_aCD_%sfilt_%s_L_%d%s_box%d_%d',con_str,filter_type,L/1000,param_num_str,box_num,year_vec(i)),'X','FFINAL','box_opt');
        fprintf('saved %d\n',year_vec(i))
        [aCd0 X{i}]
        
    end
    if plot_flag
        
        year = year_vec(i);
        filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year);

%         filename = sprintf('Qs_QL_optimization_data_%d',year_vec(i));
        load(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')
        
        lat_er = lat(box_lat);
        lon_er = lon(box_lon);
        
        t_range = 1:size(sshf_patch,3);
        
        load(sprintf('opt_aCD_%sfilt_%s_L_%d%s_box%d_%d',con_str,filter_type,L/1000,param_num_str,box_num,year_vec(i)),'X','FFINAL','box_opt');

        alpha_CD = X{i};
        
        as = alpha_CD(1);
        aL = alpha_CD(2);
        
        if numel(alpha_CD)==4
            CD_s = alpha_CD(3);
            CD_L = alpha_CD(4);
        elseif numel(alpha_CD)==3
            CD_s = alpha_CD(3);
            CD_L = alpha_CD(3);
        end
         
        sshf_model = CD_s.*(1+as.*SST_prime).*as_multiplier;
        slhf_model = CD_L.*(1+aL.*SST_prime).*aL_multiplier;
        
        sshf_model_opt = sshf_model(opt_prime_lon,opt_prime_lat,t_range);
        slhf_model_opt = slhf_model(opt_prime_lon,opt_prime_lat,t_range);
        
        ERA5_sshf = sshf_patch(opt_patch_lon,opt_patch_lat,t_range);
        ERA5_slhf = slhf_patch(opt_patch_lon,opt_patch_lat,t_range);
        
        mean_ERA5_sshf = nanmean(sshf_patch(opt_patch_lon,opt_patch_lat,t_range),3)';
        mean_ERA5_slhf = nanmean(slhf_patch(opt_patch_lon,opt_patch_lat,t_range),3)';
        
        mean_model_sshf = nanmean(sshf_model(opt_prime_lon,opt_prime_lat,t_range),3)';
        mean_model_slhf = nanmean(slhf_model(opt_prime_lon,opt_prime_lat,t_range),3)';
                
        figure(year)
        cmp_ERA5_model;
        
        cmp_ERA5_model;
        [cc(1,i),cc(2,i),rms_er(1,i),rms_er(2,i),rel_er_sshf,rel_er_slhf] = cmp_ERA5_model_corr(...
            sshf_model_opt,slhf_model_opt,ERA5_sshf,ERA5_slhf);
        
        
    end
    
    
end




if plot_flag
    figure(2)
    plot(year_vec,cc(1,:),'linewidth',2,'displayname','$$c$$ sensible')
    hold on
    plot(year_vec,cc(2,:),'linewidth',2,'displayname','$$c$$ latent')
    xlabel('year','interpreter','latex')
    ylabel('correlation coefficient','interpreter','latex')
    legend('interpreter','latex','location','best')
    set(gca,'fontsize',15)
    set(gcf,'color','w')
    update_figure_paper_size()
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_cc_%d_%s_box%d_%d',L/1000,filter_type,box_num,year),'-dpdf')
    
    figure(3)
    plot(year_vec,rms_er(1,:),'linewidth',2,'displayname','rms error sensible')
    hold on
    plot(year_vec,rms_er(2,:),'linewidth',2,'displayname','rms error latent')
    xlabel('year','interpreter','latex')
    ylabel('rms error [$$W/m^2$$]','interpreter','latex')
    legend('interpreter','latex','location','best')
    set(gca,'fontsize',15)
    set(gcf,'color','w')
    update_figure_paper_size()
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_rms_%d_%s_box%d_%d',L/1000,filter_type,box_num,year),'-dpdf')

end


