clear;clc;close all

year_vec = [2003 2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

L_vec = [500000]; % m

alpha_pos_flag = false;

%% filter set up

files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),'SST_patch','lat','lon');

m = size(files_for_size.SST_patch,1);
n = size(files_for_size.SST_patch,2);

d_lat = abs(files_for_size.lat(2)-files_for_size.lat(1));
d_lon = abs(files_for_size.lon(2)-files_for_size.lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = d_lat*m_per_deg;
dy = d_lon*m_per_deg;



% this is the same for every year
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon');

err_box_lat = [32 38];
err_box_lon = [140 160];

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));

lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

month_str = '';

for j=1:length(L_vec)
    L=L_vec(j);
    Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
    Ny = floor(L/dy)+mod(floor(L/dx),2)+1; % make Ny odd
    
    NaN_inds = isnan(files_for_size.SST_patch(:,:,1));
    
    [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
    
    for i = 1:length(year_vec)
        
        year = year_vec(i);
        dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
        load(dataFile)
        
        tt_vec = false(length(time),1);
        if ismember('D',month_str)
            tt_vec = tt_vec | month(time)==12;
        end
        if ismember('J',month_str)
            tt_vec = tt_vec | month(time)==1;
        end
        if ismember('F',month_str)
            tt_vec = tt_vec | month(time)==2;
        end
        if ismember('M',month_str)
            tt_vec = tt_vec | month(time)==3;
        end
        p = sum(tt_vec);
        tt_inds = find(tt_vec);
        
        salinity = 34*ones(m,n);% ppt for Lv calculation
                
        ZEROS = zeros(m,n,p);
        
        SST_prime = zeros(sum(patch_lon),sum(patch_lat),length(time));
        
        % coefficients
        
        load(sprintf('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha/opt_alpha_L_%d_%sCD_%d_to_%d',L/1000,con_str,2003,2007),'X');
        alpha_CD = X{year - 2002};
        
        as = alpha_CD(1);
        aL = alpha_CD(2);
        
        CD_s = alpha_CD(3);
        CD_L = alpha_CD(4);
        
        model_full_sshf = ZEROS;
        model_full_slhf = ZEROS;
        model_no_eddy_sshf = ZEROS;
        model_no_eddy_slhf = ZEROS;
        era_no_eddy_sshf = ZEROS;
        era_no_eddy_slhf = ZEROS;
        term1 = ZEROS;
        term2 = ZEROS;
        term3 = ZEROS;
        term1_CTRL = ZEROS;
        term2_CTRL = ZEROS;
        term3_CTRL = ZEROS;
        
        count = 1;
        fprintf('\n')
        for tt = tt_inds % time points
            fprintf(' processing snapshot %d of %d\n',tt,p)
            
            [SST_patch_CTRL,SST_prime] = boxcar_filter(SST_patch(:,:,tt),M);
            [P0_patch_CTRL,~] = boxcar_filter(P0_patch(:,:,tt),M);
            [q_diff_CTRL,~] = boxcar_filter(qo_patch(:,:,tt)-qa_patch(:,:,tt),M);
            
            DT_patch = SST_patch(:,:,tt) - t2m_patch(:,:,tt);
            [DT_patch_CTRL,~] = boxcar_filter(DT_patch,M);
            
            [U_mag_CTRL,~] = boxcar_filter(U_mag(:,:,tt),M);
            
            %         qo_patch_CTRL = SAM_qsatWater(SST_patch_CTRL, P0_patch(:,:,tt)) ;
            
            model_full_sshf(:,:,count) = rho_a.*c_p_air.*CD_s.*(1+as.*SST_prime).*U_mag(:,:,tt).*(DT_patch);
            model_full_slhf(:,:,count) = rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*(1+aL.*SST_prime).*U_mag(:,:,tt).*(qo_patch(:,:,tt)-qa_patch(:,:,tt));
            
            era_no_eddy_sshf(:,:,count) = boxcar_filter(sshf_patch(:,:,tt),M);
            era_no_eddy_slhf(:,:,count) = boxcar_filter(slhf_patch(:,:,tt),M);
            
            model_no_eddy_sshf(:,:,count) = rho_a.*c_p_air.*CD_s.*U_mag_CTRL.*(DT_patch_CTRL);
            model_no_eddy_slhf(:,:,count) = rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_CTRL.*(q_diff_CTRL);
            
            count = count + 1;
        end
        
        save(sprintf('model_n_ERA_data_L_%d_%s%s_%d',L/1000,con_str,month_str, year),...
            'model_full_sshf','model_full_slhf','model_no_eddy_sshf','model_no_eddy_slhf',...
            'era_no_eddy_sshf','era_no_eddy_slhf')
        
    end
end

% plot_eddy_vs_ERA;

% get_eddy_contribution_plot;
