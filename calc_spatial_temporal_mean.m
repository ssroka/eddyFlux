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
    
    load(sprintf('%sflux_terms_%d_%sfilt_%s_%d',flux_data_src,L/1000,con_str,filter_type,year))
    
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
        'lat','lon','patch_lat','patch_lon','time');
    
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
    spatial_mean = zeros(p,8,2); % sensible in the first column, latent in the second
    temporal_filter = zeros(p,8,2);
    temporal_mean = zeros(1,8,2);
    
    for j = 1:2
        
        if j == 1
            title_str = 'sensible';
        else
            title_str = 'latent';
        end
        
        for k = 1:8
            eval(sprintf('current_term = term%d;',k));
            extract_region = current_term(lon_inds,lat_inds,tt_inds,j);
            spatial_mean(:,k,j) = squeeze(nanmean(nanmean(extract_region,1),2));
            temporal_mean(1,k,j) = mean(spatial_mean(:,k,j));
            temporal_filter(:,k,j) = movmean(spatial_mean(:,k,j),temporal_filter_length_pts);
        end
        
    end
    save(sprintf('flux_terms_mean_%d_%sfilt_%s_%d_%d_%d_%d_%d',L/1000,con_str,filter_type,...
                  box_limits(1),box_limits(2),box_limits(3),box_limits(4), year),...
        'spatial_mean','temporal_mean','temporal_filter')
    fprintf('saved %d mean flux terms\n',year)
end















