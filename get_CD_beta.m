
% % % % clear;close all;clc
% % % % addpath('~/Documents/MATLAB/util/')
% % % % addpath('~/MIT/Research/eddyFlux/filter/')
% % % % addpath('~/MIT/Research/eddyFlux/ERA5_data/')
% % % % addpath('/Users/ssroka/MIT/Research/eddyFlux')
% % % % 
% % % % calc_CD_alpha_flag = true;
% % % % plot_flag = false;
% % % % 
% % % % alpha_pos_flag = false;
% % % % 
% % % % %     % Kurishio
% % % % %     lat_bnds = [25 45];
% % % % %     lon_bnds = [130 170];
% % % % 
% % % % 
% % % % L = 250000; % m
% % % % 
% % % % data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';
% % % % 
% % % % filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
% % % % 
% % % % year_vec = [2003];
% % % % 
% % % % % % 2003
% % % % % aCd0(:,1) = [  0.0105; 0.0165; 0.0014; 0.0015];
% % % % % 
% % % % % % 2004
% % % % % aCd0(:,2) = [ 0.0114; 0.0145; 0.0014; 0.0015];
% % % % % 
% % % % % % 2005
% % % % % aCd0(:,3) = [ 0.0106; 0.0126; 0.0014; 0.0015];
% % % % % 
% % % % % % 2006
% % % % % aCd0(:,4) = [ -0.0004; 0.0056; 0.0014; 0.0015];
% % % % % 
% % % % % % 2007
% % % % % aCd0(:,5) = [ 0.0125; 0.0115; 0.0014; 0.0015];
% % % % 
% % % % % aCd0 = [0.0038 0.0045 0.0014 0.0015]'; % lanczos
% % % % % aCd0 = [  0.0105; 0.0165; 0.0014; 0.0015]; % boxcar
% % % % % aCd0 = [  0.0105; 0.0165; 0.0014; 0.0015]; % boxcar
% % % % aCd0 =   [  0.0104  0.0135    0.0014]'; % fft L = 250 km 
% % % % aCd0 =   [  0.0201  0.0218    0.0014    0.0014]'; % fft L = 500 km
% % % % 
% % % % param_num_str = '_3param';
% % % % 
% % % % box_num = 1;
% % % % 
% % % % %% parameters for optimization
% % % % 
% % % % if strcmp(filter_type,'boxcar')
% % % %     box_opt = [32 38; 140 160];
% % % % elseif strcmp(filter_type,'fft')
% % % %     cf = (1/(2*L));
% % % %     %     box_opt = [32 38; 143 167];
% % % %     %     box_opt = [30 42; 144 168];
% % % % %     box_opt = [30 44.5; 148 169];
% % % % %     box_opt = [36 41.5; 143 152];
% % % % 
% % % % elseif strcmp(filter_type(1:7),'lanczos')
% % % %     cf = (1/(2*L));
% % % %     box_opt = [32 38; 140 160];
% % % % end
% % % % 
% % % % switch box_num
% % % %     case 1
% % % %         box_opt = [36 41.5; 143 152];
% % % %     case 2
% % % %         box_opt = [30 44.5; 148 169];
% % % % end
% % % % 
% % % % if alpha_pos_flag
% % % %     LB = [0;0; 0];
% % % %     UB = [Inf; Inf; Inf];
% % % %     con_str = 'cons_';
% % % % else
% % % %     LB = [-Inf; -Inf; 0];
% % % %     UB = [Inf; Inf; Inf];
% % % %     con_str = '';
% % % % end

options = optimoptions(@fmincon,'Display','iter','steptolerance',1e-5);

%% begin optimizing for CD and alpha

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'lat','lon','patch_lat','patch_lon');

file_for_box = load(sprintf('beta_Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,2003),'box_opt');

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
    
    if calc_CD_beta_flag
        
        [x,ffinal] = fmincon(@(bCd) optimize_beta_CD(bCd,year_vec(i),opt_prime_lat,opt_prime_lon,opt_patch_lat,opt_patch_lon,L,filter_type,box_num),bCd0,[],[],[],[],LB,UB,[],options);
        
        X{i} = x;
        FFINAL{i}  =ffinal;
        save(sprintf('opt_bCD_%sfilt_%s_L_%d_box%d_%d',con_str,filter_type,L/1000,box_num,year_vec(i)),'X','FFINAL','box_opt');
        fprintf('saved %d\n',year_vec(i))
        [bCd0 X{i}]
        
    end
    if plot_flag
        
        year = year_vec(i);
        filename = sprintf('beta_Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year);

%         filename = sprintf('Qs_QL_optimization_data_%d',year_vec(i));
        load(filename,'bs_multiplier','bL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch')
        
        lat_er = lat(box_lat);
        lon_er = lon(box_lon);
        
        t_range = 1:size(sshf_patch,3);
        
        load(sprintf('opt_bCD_%sfilt_%s_L_%d_box%d_%d',con_str,filter_type,L/1000,box_num,year_vec(i)),'X','FFINAL','box_opt');

        beta_CD = X{i};
        
        bs = beta_CD(1);
        bL = beta_CD(2);
        
        if numel(beta_CD)==4
            CD_s = beta_CD(3);
            CD_L = beta_CD(4);
        elseif numel(beta_CD)==3
            CD_s = beta_CD(3);
            CD_L = beta_CD(3);
        end
         
        sshf_model = CD_s.*(U_bar+bs.*SST_prime).*bs_multiplier;
        slhf_model = CD_L.*(U_bar+bL.*SST_prime).*bL_multiplier;
        
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


