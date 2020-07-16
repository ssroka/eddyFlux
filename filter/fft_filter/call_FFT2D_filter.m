% =========================================================================
% call_FFT2D_filter
% =========================================================================
% Created By    : Sydney Sroka
% Last Edited By: Sydney Sroka
% Date          : 07/14/2020
%
% ------------------------
% This script was written to practice calling FFT2D_filter.

clear;close all;clc

cT = 1000000; % cutoff period;

year = [2003];

box_limits = [30 42; 148 168];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

%% begin

% ------- load and set up 2D data -----------------------------------
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year));

lat = double(lat);
lon = double(lon);

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

lat_inds = lat_er>=box_limits(1,1) & lat_er<=box_limits(1,2);
lon_inds = lon_er>=box_limits(2,1) & lon_er<=box_limits(2,2);

lat_er = lat_er(lat_inds);
lon_er = lon_er(lon_inds);

s = nanmean(SST_patch(lon_inds,lat_inds,:),3)';

m = size(SST_patch,1);
n = size(SST_patch,2);

d_lat = abs(lat(2)-lat(1));
d_lon = abs(lon(2)-lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

cf = 1/cT;

figure
contourf(lon_er,lat_er,nanmean(slhf_patch(lon_inds,lat_inds,:),3)')

% ------- preprocess data for filtering -----------------------------------
s = s - mean(s(:));

s = detrend(s); % detrend columns
s = detrend(s')'; % detrend rows

debug_flag = true;

[S_LP,S_HP] = FFT2D_filter(s,dx,cf,debug_flag,lon_er,lat_er);

set(gcf,'position',[440           1        1001         797])
update_figure_paper_size()
print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/test_fft_%d',year),'-dpdf')


