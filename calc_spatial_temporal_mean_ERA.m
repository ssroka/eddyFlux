clear;clc;close all

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';
flux_data_src = '/Volumes/SydneySroka_Remy/eddyFlux_data/old_lanczos_filtered_results/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/lanczosfilter/')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

L = 500000; % m

intvl = 1; % look at every intvl'th timpepoint

plt_error_box = false;

year_vec = [2003:2019];

add_ssh_flag = false;
%
month_str = 'DJFM';

alpha_pos_flag = false;

filter_type = 'lanczos'; % filter type 'lanczos' or 'boxcar'

box_limits = [34 41; 143 168];

temporal_filter_length = 240;% hours (equiv to 10 days), min is 6

%% begin

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end
for i = 1:length(year_vec)
    
    year = year_vec(i);
    
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
        'lat','lon','patch_lat','patch_lon','time',...
        'slhf_patch','sshf_patch');
    
    temporal_filter_length_pts = ceil(temporal_filter_length/hours(time(2)-time(1)));
    
    lat_er = lat(patch_lat);
    lon_er = lon(patch_lon);
    
    lat_inds = lat_er>=box_limits(1,1) & lat_er<=box_limits(1,2);
    lon_inds = lon_er>=box_limits(2,1) & lon_er<=box_limits(2,2);
    
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
    spatial_mean_L = zeros(p,1); % sensible in the first column, latent in the second
    spatial_mean_S = zeros(p,1); % sensible in the first column, latent in the second
    temporal_filter_L = zeros(p,1);
    temporal_filter_S = zeros(p,1);
    temporal_mean_L = zeros(1,1);
    temporal_mean_S = zeros(1,1);
    
    
    extract_latent = slhf_patch(lon_inds,lat_inds,tt_inds,:);
    extract_sensible = sshf_patch(lon_inds,lat_inds,tt_inds,:);
    
    spatial_mean_L(:,1) = squeeze(nanmean(nanmean(extract_latent,1),2));
    spatial_mean_S(:,1) = squeeze(nanmean(nanmean(extract_sensible,1),2));
    
    temporal_mean_L(1,1) = mean(spatial_mean_L(:,1));
    temporal_mean_S(1,1) = mean(spatial_mean_S(:,1));
    
    temporal_filter_L(:,1) = movmean(spatial_mean_L(:,1),temporal_filter_length_pts);
    temporal_filter_S(:,1) = movmean(spatial_mean_S(:,1),temporal_filter_length_pts);
    
    save(sprintf('flux_terms_mean_ERA_%d_%d_%d_%d_%d',...
        box_limits(1),box_limits(3),box_limits(2),box_limits(4), year),...
        'spatial_mean_L','spatial_mean_S','temporal_mean_L','temporal_mean_S',...
        'temporal_filter_L','temporal_filter_S')
    
    fprintf('saved %d mean ERA flux terms\n',year)
end















