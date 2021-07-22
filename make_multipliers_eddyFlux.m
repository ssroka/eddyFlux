function [] = make_multipliers_eddyFlux(year_vec,L,data_src,filter_type,box_num,box_opt,cf,dx,debug_flag,model_str)
%% begin optimizing for CD and alpha
files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'SST_patch','lat','lon','patch_lat','patch_lon');

setup_lat_lon_vec;
constant_vals = load('env_const');
for i = 1:length(year_vec)
    year = year_vec(i);
    clearvars -except year M data_src m n constant_vals L filter_type box_opt...
        dx cf year_vec lat_patch_2_box_TF lon_patch_2_box_TF lat_box lon_box debug_flag box_num model_str
    
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile,'time','SST_patch','U_mag','t2m_patch',...
        'qo_patch','qa_patch','sshf_patch','slhf_patch')
            
    p = size(SST_patch,3); % needs to be recalculated because of leap year
    salinity = 34*ones(m,n,p);% ppt for Lv calculation
    
    SST_prime = zeros(m,n,p);
    U_bar = zeros(m,n,p);
    
    for tt = 1:length(time)
            [SST_bar] = FFT2D_filter(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)',dx,cf,debug_flag,lon_box,lat_box);
            SST_prime(:,:,tt) = SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt) - SST_bar';
            [U_bar_tt] = FFT2D_filter(U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)',dx,cf,debug_flag,lon_box,lat_box);
            U_bar(:,:,tt) = U_bar_tt';
            fprintf('getting SST prime and U_bar for time step %d of %d\n',tt,p)
    end
    
    filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year);
    switch model_str
        case 'alpha'
        as_multiplier = constant_vals.rho_a.*constant_vals.c_p_air.*U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,:).*(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)-t2m_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:));
        aL_multiplier = constant_vals.rho_a.*SW_LatentHeat(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),'K',salinity,'ppt').*U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,:).*(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:));
        save(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch','L','box_opt')
        
        case 'beta'
        bs_multiplier = constant_vals.rho_a.*constant_vals.c_p_air.*(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)-t2m_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:));
        bL_multiplier = constant_vals.rho_a.*SW_LatentHeat(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),'K',salinity,'ppt').*(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:));
        save(filename,'bs_multiplier','bL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch','L','box_opt')
        
        case 'alphabeta'
        abs_multiplier = constant_vals.rho_a.*constant_vals.c_p_air.*(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)-t2m_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:));
        abL_multiplier = constant_vals.rho_a.*SW_LatentHeat(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),'K',salinity,'ppt').*(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:));
        save(filename,'abs_multiplier','abL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch','L','box_opt')
        
    end
    
    fprintf('saved %d\n',year)
end
