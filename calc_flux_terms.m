clear;clc;close all

year_vec = [2003 2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')

L = 500000; % m

intvl = 1; % look at every intvl'th timpepoint

alpha_pos_flag = false;

%% filter set up

files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),'SST_patch','lat','lon','patch_lat','patch_lon');

m = size(files_for_size.SST_patch,1);
n = size(files_for_size.SST_patch,2);

d_lat = abs(files_for_size.lat(2)-files_for_size.lat(1));
d_lon = abs(files_for_size.lon(2)-files_for_size.lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = d_lat*m_per_deg;
dy = d_lon*m_per_deg;

Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
Ny = floor(L/dy)+mod(floor(L/dx),2)+1; % make Ny odd

NaN_inds = isnan(files_for_size.SST_patch(:,:,1));

[M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);

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

salinity = 34*ones(m,n);% ppt for Lv calculation

for i = 1:length(year_vec)
    year = year_vec(i);
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile)
    
    SST_prime = zeros(sum(patch_lon),sum(patch_lat),length(time));
    
    % coefficients
    load(sprintf('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha/opt_alpha_L_%d_%sCD_%d_to_%d',L/1000,con_str,2003,2007),'X');
    alpha_CD = X{year - 2002};
    
    as = alpha_CD(1);
    aL = alpha_CD(2);
    
    CD_s = alpha_CD(3);
    CD_L = alpha_CD(4);
    
    p = size(files_for_size.SST_patch,3); % re calculate for leap year
        
    ZEROS = zeros(m,n,p,2);
    
    
    term1 = ZEROS;
    term2 = ZEROS;
    term3 = ZEROS;
    term4 = ZEROS;
    term5 = ZEROS;
    term6 = ZEROS;
    term7 = ZEROS;
    term8 = ZEROS;

    
    count = 1;
    fprintf('\n')
    for tt = 1:intvl:p % time points
        fprintf(' processing snapshot %d of %d\n',tt,p)
        
        [SST_patch_CTRL,SST_prime] = boxcar_filter(SST_patch(:,:,tt),M);
        [P0_patch_CTRL,~] = boxcar_filter(P0_patch(:,:,tt),M);
        
        [q_diff_CTRL,~] = boxcar_filter(qo_patch(:,:,tt)-qa_patch(:,:,tt),M);
        q_diff_prime = (qo_patch(:,:,tt)-qa_patch(:,:,tt)) - q_diff_CTRL;

        DT_patch = SST_patch(:,:,tt) - t2m_patch(:,:,tt);
        [DT_patch_CTRL,~] = boxcar_filter(DT_patch,M);
        DT_diff_prime = (DT_patch) - DT_patch_CTRL;

        
        [U_mag_CTRL,~] = boxcar_filter(U_mag(:,:,tt),M);
        U_mag_prime = U_mag(:,:,tt) - U_mag_CTRL;
        
        
        term1(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*U_mag_CTRL.*DT_patch_CTRL;
        
        term2(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*as.*DT_patch_CTRL.*SST_prime.*U_mag_prime;
        term3(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*as.*SST_prime.*U_mag_CTRL.*DT_diff_prime;
        term4(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*U_mag_prime.*DT_diff_prime;
        term5(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*as.*SST_prime.*U_mag_prime.*DT_diff_prime;
        
        term6(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*as.*SST_prime.*U_mag_CTRL.*DT_patch_CTRL;
        term7(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*U_mag_prime.*DT_patch_CTRL;
        term8(:,:,count,1) =  rho_a.*c_p_air.*CD_s.*U_mag_CTRL.*DT_diff_prime;
                
        term1(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_CTRL.*q_diff_CTRL;
        
        term2(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*aL.*q_diff_CTRL.*SST_prime.*U_mag_prime;
        term3(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*aL.*SST_prime.*U_mag_CTRL.*q_diff_prime;
        term4(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_prime.*q_diff_prime;
        term5(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*aL.*SST_prime.*U_mag_prime.*q_diff_prime;
        
        term6(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*aL.*SST_prime.*U_mag_CTRL.*q_diff_CTRL;
        term7(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_prime.*q_diff_CTRL;
        term8(:,:,count,2) =  rho_a.*SW_LatentHeat(SST_patch(:,:,tt),'K',salinity,'ppt').*CD_L.*U_mag_CTRL.*q_diff_prime;
                
        count = count + 1;
    end
        save(sprintf('flux_terms_%d_%s%d',L/1000,con_str,year),...
        'term1','term2','term3','term4','term5','term6','term7','term8')
    
end

;
% get_eddy_contribution_plot;
