clc;clear;close all
addpath('~/Documents/MATLAB/util/')

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'point'; % 'box_sum' 'box_mean' or 'point'
err_box_lat = [32 38];
err_box_lon = [140 160];

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

L = 500000;%m

% as_multiplier = rho_a.*c_p_air.*U_mag.*DT_patch;
% aL_multiplier = rho_a.*Lv.*U_mag.*(qo_patch-qa_patch);

for year = 2004
    
    
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(filename_data,'lat','lon','patch_lat','patch_lon','rho_a','c_p_air',...
        'U_mag','DT_patch','Lv','qo_patch','qa_patch','SST_prime','sshf_patch',...
        'slhf_patch','time','SST_patch','P0_patch','t2m_patch','RH_patch')
    
    patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));
    
    err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
    err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
    
    %     6 params
    %     load X_FFINAL_2007.mat
%     aCd0 = X{year-2002};
%     
%     CD_ref_s = aCd0(5);
%     CD_ref_L = aCd0(6);
%     
%     as_est =@(SST_prime) aCd0(1)+aCd0(2)*abs(SST_prime);
%     aL_est =@(SST_prime) aCd0(3)+aCd0(4)*abs(SST_prime);
    
%   4 params
    load(sprintf('X_FFINAL_%d_const_only',year),'X');
    aCd0 = X{year-2002};
    CD_ref_s = aCd0(3);
    CD_ref_L = aCd0(4);
    
    as_est =@(SST_prime) aCd0(1);
    aL_est =@(SST_prime) aCd0(2);
    
    %     filename_mult = sprintf('mult_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year-2000,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    %     load(filename_mult,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')
    
%     sshf_as_constant_full = CD_ref_s.*(1+as_est(SST_prime).*SST_prime).*as_multiplier;
%     slhf_aL_constant_full = CD_ref_L.*(1+aL_est(SST_prime).*SST_prime).*aL_multiplier;
%     
%     
%     sshf_as_constant_eddy = sshf_as_constant_full - sshf_as_constant_no_eddy;
%     slhf_aL_constant_eddy = slhf_aL_constant_full - slhf_aL_constant_no_eddy;
%     
%     tot_flux_1 = sshf_as_constant + slhf_aL_constant;
    
    sshf_eddy = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_eddy = zeros(sum(patch_lon),sum(patch_lat),length(time));
    sshf_sm = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_sm = zeros(sum(patch_lon),sum(patch_lat),length(time));
    sshf_diff = zeros(sum(patch_lon),sum(patch_lat),length(time));
    slhf_diff = zeros(sum(patch_lon),sum(patch_lat),length(time));
%     SST_smooth =  zeros(sum(patch_lon),sum(patch_lat),length(time));
%     SST_prime_smooth  =  zeros(sum(patch_lon),sum(patch_lat),length(time));
    term1  =  zeros(sum(patch_lon),sum(patch_lat),length(time));
    term2  =  zeros(sum(patch_lon),sum(patch_lat),length(time));
    term3  =  zeros(sum(patch_lon),sum(patch_lat),length(time));
    
%     [X,Y] = create_grid(lon(patch_lon),lat(patch_lat));
%     X_dist = abs(X(end,end))-(X(end,1));
%     Y_dist = abs(Y(1,1)-Y(end,1));
%     Nx     = floor(X_dist/L);
%     Ny     = floor(Y_dist/L);
    
    Nx = 19;
    Ny = 19;

[M] = my_smoother_mat(sum(patch_lon),sum(patch_lat),Ny,Nx);
    
    for tt = 1:length(time)
        
        % ---- qo ----
        qo_patch = SAM_qsatWater(SST_patch(:,:,tt), P0_patch(:,:,tt)) ;
        
        % ---- qa ----
        e_sat = SAM_psatWater(t2m_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat;
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch = r./(1+r);
        
        dq_eddy = qo_patch - qa_patch;
        dq_no_eddy = smooth_mat(qo_patch - qa_patch,M);
        dq_no_eddy2 = my_smoother(qo_patch - qa_patch,Ny,Nx);
        Umag_no_eddy = smooth_mat(U_mag(:,:,tt),M);
        
        % add the smooth U_mag 
        
        DT_eddy = SST_patch(:,:,tt)-t2m_patch(:,:,tt);
        DT_no_eddy = smooth_mat(SST_patch(:,:,tt)-t2m_patch(:,:,tt),M);
        T_o_prime = SST_patch(:,:,tt) - smooth_mat(SST_patch(:,:,tt),M);
        T_o_prime_sq_sm = smooth_mat(T_o_prime.^2,M);
        
        T_a_prime = t2m_patch(:,:,tt) - smooth_mat(t2m_patch(:,:,tt),M);
        T_a_prime_T_o_prime_sm = smooth_mat(T_a_prime.*T_o_prime,M);
        
        sshf_eddy(:,:,tt) = rho_a.*c_p_air.*CD_ref_s.*(1+as_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*(DT_eddy);
        slhf_eddy(:,:,tt) = rho_a.*Lv.*CD_ref_L.*(1+aL_est(SST_prime(:,:,tt)).*SST_prime(:,:,tt)).*U_mag(:,:,tt).*(dq_eddy);
        
        sshf_sm(:,:,tt) = rho_a.*c_p_air.*CD_ref_s.*Umag_no_eddy.*(DT_no_eddy);
        slhf_sm(:,:,tt) = rho_a.*Lv.*CD_ref_L.*Umag_no_eddy.*(dq_no_eddy);
        
        sshf_diff(:,:,tt) = smooth_mat(sshf_eddy(:,:,tt) - sshf_sm(:,:,tt),M);
        slhf_diff(:,:,tt) = smooth_mat(slhf_eddy(:,:,tt) - slhf_sm(:,:,tt),M);
        
        term1(:,:,tt) = rho_a.*c_p_air.*CD_ref_s.*Umag_no_eddy.*DT_no_eddy;
        term2(:,:,tt) = rho_a.*c_p_air.*CD_ref_s.*as_est(SST_prime(:,:,tt)).*Umag_no_eddy.*T_o_prime_sq_sm;
        term3(:,:,tt) = -rho_a.*c_p_air.*CD_ref_s.*as_est(SST_prime(:,:,tt)).*Umag_no_eddy.*T_a_prime_T_o_prime_sm;
        
    end
    close all
    plot_model_eddy_min_no_eddy
    plot_terms_fm_deriv
    save(sprintf('%d_model_eddy_min_no_eddy_results.mat',year))
end















