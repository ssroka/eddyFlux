clear;clc;close all

year_vec = [2003];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

L_vec = [500000]; % m

alpha_pos_flag = false;

month_str = 'DJFM';

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
                
        SST_prime = zeros(sum(patch_lon),sum(patch_lat),length(time));
        
        fprintf('\n')
        for tt = tt_inds(1:50:end)' % time points
            fprintf(' processing snapshot %d of %d\n',tt,tt_inds(end))
            
            [q_diff_CTRL,q_diff_prime] = boxcar_filter(qo_patch(:,:,tt)-qa_patch(:,:,tt),M);
            [q_o_CTRL,q_o_prime] = boxcar_filter(qo_patch(:,:,tt),M);
            [q_a_CTRL,q_a_prime] = boxcar_filter(qa_patch(:,:,tt),M);
            
            [DT_patch_CTRL,DT_patch_prime] = boxcar_filter(DT_patch(:,:,tt),M);
            
            figure(1)
            plot(q_diff_prime-(q_o_prime-q_a_prime),'ko')
            title('$\Delta q'' - (q_o''-q_a'')','interpreter','latex')
            hold on
            figure(2)
            plot(q_diff_prime,DT_patch_prime,'k*')
            xlabel('$\Delta q''$','interpreter','latex')
            ylabel('$\Delta T''$','interpreter','latex')
            hold on
            figure(3)
            plot(qo_patch(:,:,tt)-qa_patch(:,:,tt),DT_patch(:,:,tt),'k*')
            hold on
            xlabel('(q_o-q_a)','interpreter','latex')
            ylabel('(T_o-T_a)','interpreter','latex')
            figure(4)
            plot(q_diff_CTRL,DT_patch_CTRL,'k*')
            hold on
            xlabel('$\overline{(q_o-q_a)}$','interpreter','latex')
            ylabel('$\overline{(T_o-T_a)}$','interpreter','latex')
        end

        
    end
end

% plot_eddy_vs_ERA;

% get_eddy_contribution_plot;
