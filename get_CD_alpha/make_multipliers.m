% as_multiplier = rho_a.*c_p_air.*U_mag.*(SST_patch-t2m_patch);
% aL_multiplier = rho_a.*Lv.*U_mag.*(qo_patch-qa_patch);


clear;close all;clc
addpath('~/Documents/MATLAB/util/')
addpath('~/MIT/Research/eddyFlux/filter/')
addpath('~/MIT/Research/eddyFlux/ERA5_data/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

%     % Kurishio
%     lat_bnds = [25 45];
%     lon_bnds = [130 170];


data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

L = 800000; % m

%% begin optimizing for CD and alpha
files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),'SST_patch','lat','lon');

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


for year = 2003:2007
    
    clearvars -except year M data_src L m n
    
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile)
    
    p = size(SST_patch,3); % needs to be recalculated because of leap year
    salinity = 34*ones(m,n,p);% ppt for Lv calculation

    SST_prime = zeros(size(SST_patch,1),size(SST_patch,2),size(SST_patch,3));
    
    for tt = 1:length(time)
        [~,SST_prime(:,:,tt)] = boxcar_filter(SST_patch(:,:,tt),M);
        fprintf('getting SST prime for time step %d of %d\n',tt,p)
    end
    
    as_multiplier = rho_a.*c_p_air.*U_mag.*(SST_patch-t2m_patch);
    aL_multiplier = rho_a.*SW_LatentHeat(SST_patch,'K',salinity,'ppt').*U_mag.*(qo_patch-qa_patch);
    
    filename = sprintf('Qs_QL_optimization_data_L_%d_%d',L/1000,year);
    
    save(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch','L')
    fprintf('saved %d\n',year)
end





























































