ccc

reanalysis_src = 'NCEP';
year = 2003;

if strcmp(reanalysis_src,'NCEP')
lat_name = 'LAT461_54033_68';
lon_name = 'LON541_70029_120';
else
lat_name = 'LATITUDE461_54133_69';
lon_name = 'LONGITUDE541_70129_121';
end



cd(sprintf('~/Downloads/%s',reanalysis_src))
figure(1)
subplot(1,2,1)
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
debug_flag = false;

lon = ncread(sprintf('SSTA_%d.nc',year),lon_name);
lat = ncread(sprintf('SSTA_%d.nc',year),lat_name);
SST = ncread(sprintf('SSTA_%d.nc',year),'ANOM');

m = size(SST,1);
n = size(SST,2);
p = size(SST,3);

% only used to get the fft filter sampling rate;
d_lat = abs(lat(2)-lat(1));
m_per_deg = 111320;
dx = abs(d_lat*m_per_deg);

if strcmp(filter_type,'fft')
    cf = (1/(2*L));
end

SST_prime = zeros(m,n,p);
nan_matrix = isnan(SST);
SST(find(nan_matrix)) = 0.5;
for tt = 1:p
[SST_bar(:,:,tt)] = FFT2D_filter(SST(:,:,tt)',dx,cf,debug_flag,lon,lat);
SST_prime(:,:,tt) = SST(:,:,tt) - SST_bar(:,:,tt)';
end
imagesc(lon,lat, mean(SST_prime,3)')
title(sprintf('%d %s FFT filter',year,reanalysis_src))
colorbar
set(gca,'clim',[-2 2])

%% plot HPF
subplot(1,2,2)
lon = ncread(sprintf('SSTHP_%d.nc',year),lon_name);
lat = ncread(sprintf('SSTHP_%d.nc',year),lat_name);
SST_HP = ncread(sprintf('SSTHP_%d.nc',year),'SST_HP');
imagesc(lon,lat, mean(SST_HP,3)')
title(sprintf('%d %s High pass boxcar',year,reanalysis_src))
colorbar
set(gca,'clim',[-2 2])
set(gcf,'position',[1         433        1416         365])
update_figure_paper_size()
print(sprintf('%s_%d_fm_Soumi',reanalysis_src,year),'-dpdf')

