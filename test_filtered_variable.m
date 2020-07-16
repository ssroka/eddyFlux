clear;clc;close all

year_vec = [2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')
addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/lanczosfilter')

L_vec = [200000]; % m

alpha_pos_flag = false;

month_str = 'DJFM';

filter_type = 'lanczos_xy'; % filter type 'lanczos' or 'boxcar'


%% filter set up
load('env_const.mat')

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'SST_patch','lat','lon','patch_lat','patch_lon');

m = size(SST_patch,1);
n = size(SST_patch,2);

d_lat = abs(lat(2)-lat(1));
d_lon = abs(lon(2)-lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

% this is the same for every year
% load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
%     'lat','lon','patch_lat','patch_lon');

err_box_lat = [32 38];
err_box_lon = [140 160];

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));

lat_er = lat(patch_lat);
% lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
% lon_er = lon_er(err_box_bnds_lon);

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end


for j=1:length(L_vec)
    L=L_vec(j);
    Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
    Ny = floor(L/dy)+mod(floor(L/dx),2)+1; % make Ny odd
    
    NaN_inds = isnan(SST_patch(:,:,1));
    
    if strcmp(filter_type,'boxcar')
        [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
    elseif strcmp(filter_type(1:7),'lanczos')
        cf = (1/(2*L));
    end
    
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
        
        
        count = 1;
        fprintf('\n')
%         for tt = tt_inds' % time points
%             fprintf(' processing snapshot %d of %d\n',tt,tt_inds(end))
            
            if strcmp(filter_type,'boxcar')
                [SST_patch_CTRL,SST_prime] = boxcar_filter(SST_patch(:,:,tt),M);
                %             [P0_patch_CTRL,~] = boxcar_filter(P0_patch(:,:,tt),M);
                [q_diff_CTRL,~] = boxcar_filter(qo_patch(:,:,tt)-qa_patch(:,:,tt),M);
                [DT_patch_CTRL,~] = boxcar_filter(DT_patch(:,:,tt),M);
                [U_mag_CTRL,~] = boxcar_filter(U_mag(:,:,tt),M);
                era_no_eddy_sshf(:,:,count) = boxcar_filter(sshf_patch(:,:,tt),M);
                era_no_eddy_slhf(:,:,count) = boxcar_filter(slhf_patch(:,:,tt),M);
            elseif strcmp(filter_type(1:7),'lanczos')
                A = nanmean(SST_patch(:,:,:),3);
                [A_sm,A_prime,A_sm_x,A_sm_y]= lanczos_filter(A,dx,cf);
            end
            
            plot_filtered_variable;
%         end

        
    end
end

% plot_eddy_vs_ERA;

% get_eddy_contribution_plot;
