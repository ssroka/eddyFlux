% as_multiplier = rho_a.*c_p_air.*U_mag.*(SST_patch-t2m_patch);
% aL_multiplier = rho_a.*Lv.*U_mag.*(qo_patch-qa_patch);


clear;close all;clc
addpath('~/Documents/MATLAB/util/')
addpath('~/MIT/Research/eddyFlux/filter/')
addpath('~/MIT/Research/eddyFlux/filter/lanczosfilter/')
addpath('~/MIT/Research/eddyFlux/filter/fft_filter/')
addpath('~/MIT/Research/eddyFlux/')
addpath('~/MIT/Research/eddyFlux/ERA5_data/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

%     % Kurishio
%     lat_bnds = [25 45];
%     lon_bnds = [130 170];

year_vec = [2003 2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

L = 500000; % m

filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'

%% begin optimizing for CD and alpha
files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'SST_patch','lat','lon','patch_lat','patch_lon');

if strcmp(filter_type,'boxcar')
    box_limits = [25 45; 130 170];
elseif strcmp(filter_type,'fft')
    cf = (1/(2*L));
    box_limits = [30 42; 148 168];
elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
    box_limits = [25 45; 130 170];
end

box_lat = files_for_size.lat>=box_limits(1,1) & files_for_size.lat<=box_limits(1,2);
box_lon = files_for_size.lon>=box_limits(2,1) & files_for_size.lon<=box_limits(2,2);

lat_er = files_for_size.lat(box_lat);
lon_er = files_for_size.lon(box_lon);

% to index out of *_patch fields
inds_lat = box_lat(files_for_size.patch_lat);
inds_lon = box_lon(files_for_size.patch_lon);

n = length(lat_er);
m = length(lon_er);

d_lat = abs(files_for_size.lat(2)-files_for_size.lat(1));
d_lon = abs(files_for_size.lon(2)-files_for_size.lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
Ny = floor(L/dy)+mod(floor(L/dy),2)+1; % make Ny odd

NaN_inds = isnan(files_for_size.SST_patch(:,:,1));

load('env_const.mat'); % load rho_a and c_p_air

if strcmp(filter_type,'boxcar')
    [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
end

for i = 1:length(year_vec)
    year = year_vec(i);
    clearvars -except year M data_src m n rho_a L c_p_air filter_type...
        dx cf year_vec box_limits inds_lat inds_lon
    
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile)
            
    p = size(SST_patch,3); % needs to be recalculated because of leap year
    salinity = 34*ones(m,n,p);% ppt for Lv calculation
    
    SST_prime = zeros(m,n,p);
    
    for tt = 1:length(time)
        if strcmp(filter_type,'boxcar')
            [~,SST_prime(:,:,tt)] = boxcar_filter(SST_patch(inds_lon,inds_lat,tt),M);
        elseif strcmp(filter_type,'fft')
            [SST_patch_CTRL] = FFT2D_filter(SST_patch(inds_lon,inds_lat,tt),dx,cf);
            SST_prime(:,:,tt) = SST_patch(inds_lon,inds_lat,tt) - SST_patch_CTRL;
        elseif strcmp(filter_type(1:7),'lanczos')
            [~,SST_prime(:,:,tt)] = lanczos_filter(SST_patch(inds_lon,inds_lat,tt),dx,cf);
        end
        fprintf('getting SST prime for time step %d of %d\n',tt,p)
    end
    
    as_multiplier = rho_a.*c_p_air.*U_mag(inds_lon,inds_lat,:).*(SST_patch(inds_lon,inds_lat,:)-t2m_patch(inds_lon,inds_lat,:));
    aL_multiplier = rho_a.*SW_LatentHeat(SST_patch(inds_lon,inds_lat,:),'K',salinity,'ppt').*U_mag(inds_lon,inds_lat,:).*(qo_patch(inds_lon,inds_lat,:)-qa_patch(inds_lon,inds_lat,:));
    
    filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_%d',L/1000,filter_type,year);
    
    save(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch','L','box_limits')
    fprintf('saved %d\n',year)
end



