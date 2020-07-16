
ccc

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/')

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

year = 2003;

L = 500000;

%% begin

dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
load(dataFile)


lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

d_lat = (lat_er(2)-lat_er(1));
d_lon = (lon_er(2)-lon_er(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

% L_deg = L/m_per_deg;
% cf = 1/(2*L_deg);

subplot(2,3,1)
[~,h] = contourf(lon_er,lat_er,nanmean(SST_patch,3)');
title(sprintf('ERA5 full SST $$K$$'),'interpreter','latex')

[SST_LP,SST_prime] = fft_filter(nanmean(SST_patch,3)',dx,cf);

subplot(2,3,2)
[~,h] = contourf(lon_er,lat_er,SST_LP);
title(sprintf('Lanczos filter L = 500km (2L cutoff)'),'interpreter','latex')

subplot(2,3,3)
[~,h] = contourf(lon_er,lat_er,SST_prime);
title(sprintf('full minus Lanczos'),'interpreter','latex')

m = size(SST_patch,1);
n = size(SST_patch,2);

Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
Ny = floor(L/dy)+mod(floor(L/dy),2)+1; % make Ny odd

NaN_inds = isnan(SST_patch(:,:,1));

[M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
subplot(2,3,4)
[~,h] = contourf(lon_er,lat_er,nanmean(SST_patch,3)');
title(sprintf('ERA5 full SST $$K$$'),'interpreter','latex')

[SST_patch_CNTL,SST_patch_prime] = boxcar_filter(nanmean(SST_patch,3),M);


subplot(2,3,5)
[~,h] = contourf(lon_er,lat_er,SST_patch_CNTL');
title(sprintf('boxcar filter L = 500km'),'interpreter','latex')

subplot(2,3,6)
[~,h] = contourf(lon_er,lat_er,SST_patch_prime');
title(sprintf('full minus boxcar filter L = 500km '),'interpreter','latex')


for i = 1:6
    subplot(2,3,i)
    set(gca,'fontsize',25)
end
set(gcf,'position',[25           1        1416         787],'color','w')


figure(1)
update_figure_paper_size()
print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/test_fft_%d_%d',L/1000,year),'-dpdf')





















