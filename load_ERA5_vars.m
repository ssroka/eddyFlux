clear
% load('/Volumes/SydneySroka_Anton/winter_02_lat_lon_time.mat')
% load('/Volumes/SydneySroka_Anton/winter_02_heatFlux.mat')
% load('/Volumes/SydneySroka_Anton/winter_15_lat_lon_time.mat')
% load('/Volumes/SydneySroka_Anton/winter_15_heatFlux.mat')
% load('/Volumes/SydneySroka_Anton/winter_02_v.mat')

src = '/Volumes/SydneySroka_Anton/adaptor.mars.internal-1562361444.2138896-21299-5-befc5b3b-babe-4a4d-ab28-cdccbdc58cb8.nc';

% SRC_INFO = ncinfo(src);
% SRC_INFO.Variables.Name

time = ncread(src,'time'); % hours


topLeft = [40,(360 -70)];
botRight = [30,(360 -25)];

lat_tot = ncread(src,'latitude');
lon_tot = ncread(src,'longitude');

lat_inds = (lat_tot<=topLeft(1) & lat_tot>=botRight(1));
lon_inds = (lon_tot>=topLeft(2) & lon_tot<=botRight(2));

lat = lat_tot(lat_inds);
lon = lon_tot(lon_inds);

% convert lat/lon into meters
R = 6371000;% radius in meters

y_km = linspace(0,1223,length(lat));
x_km = linspace(0,3261,length(lon));

U10 = ncread(src,'u10');
U10 = mean(U10(lon_inds,lat_inds,:),3)';

SST = ncread(src,'sst');
SST_field = mean(SST(lon_inds,lat_inds,1:4*7*4),3)';
imagesc(x_km,y_km,SST_field);
xlabel('km')
ylabel('km')

figure 
imagesc(lon_tot,lat_tot,mean(SST(:,:,1:4:7),3)')
set(gca,'ydir','normal')
xlabel('km')
ylabel('km')

V10 = ncread(src,'v10');
V10 = mean(V10(lon_inds,lat_inds,:),3)';

P = ncread(src,'sp');
P = mean(P(lon_inds,lat_inds,:),3)';


SST_CTRL = imgaussfilt(SST_field, 10);
close all
t2m_lpf = imgaussfilt(mean(t2m,3), 10);
imagesc(mean(t2m,3)'-t2m_lpf');colorbar
title('air')
figure
sst_lpf = imgaussfilt(sst, 10);
imagesc(sst'-sst_lpf');colorbar

hold on
figure
imagesc(sshf')

T_prime = SST_field - SST_CTRL;

DT = SST_field-mean(t2m,3)';

U = mean(u10,3)';
V = mean(v10,3)';
P = mean(sp,3)';