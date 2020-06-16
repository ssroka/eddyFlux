clear;close all;clc

plt_error_box = false;

t_range = 1:484;

year = 2007;

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

add_ssh_flag = false;

L = 500000 ; % m

alpha_pos_flag = true;

%% begin
if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year));

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

m = size(SST_patch,1);
n = size(SST_patch,2);

d_lat = abs(lat(2)-lat(1));
d_lon = abs(lon(2)-lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = d_lat*m_per_deg;
dy = d_lon*m_per_deg;

Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
Ny = floor(L/dy)+mod(floor(L/dx),2)+1; % make Ny odd

NaN_inds = isnan(SST_patch(:,:,1));

[M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);

[SST_patch_CTRL,SST_prime] = boxcar_filter(SST_patch(:,:,1),M);

fft_region = [143 169; 29.5 41];

fft_patch_lat = lat>fft_region(1,1) & lat<fft_region(1,2) ;
fft_patch_lon = lon>fft_region(2,1) & lon<fft_region(2,2) ;



ax = figure(1);
[~,h] = contourf(lon_er,lat_er,SST_prime');
title(sprintf('SST $$[K]$$'),'interpreter','latex')
format_fig(h,ax)


function [] = format_fig(h,plt_num,max_val,min_val)

set(h,'edgecolor','none')
set(gca,'ydir','normal','fontsize',15)
colorbar
xlabel('deg')
ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end