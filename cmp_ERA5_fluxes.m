clear;clc;close all
% This file is to look at the sensible and latent heat variation between
% different winters in ERA 5 data
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';
% src = '/Volumes/SydneySroka_Anton/ERA5_2017_OND.nc';
% src = '/Users/ssroka/Downloads/adaptor.mars.internal-1560194023.5105438-13354-5-264c1c9e-f2f3-4535-92cc-3a69e7482cb8.nc';

for year = 1:2
    
    switch year
        case 1 % 2003
            
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
        case 2 % 2007
            
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2007_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2006_Dec_31d_6h.nc';
    end

time = double(ncread(srcJFM,'time'))*60*60; % hours
time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');

lsm = ncread(land_src,'lsm');

lat = ncread(srcJFM,'latitude');
lon = ncread(srcJFM,'longitude');

[X,Y] = create_grid(lat,lon);

slhf = nanmean(cat(3,ncread(srcD,'slhf'),ncread(srcJFM,'slhf'))./(60*60),3); % J m^-2 / (24 hours in seconds) = W m^-2
sshf = nanmean(cat(3,ncread(srcD,'sshf'),ncread(srcJFM,'sshf'))./(60*60),3); % J m^-2 / (24 hours in seconds) = W m^-2

% plot a snapshot of latent heat flux

GS_lat = (lat>25)&(lat<45);
GS_lon = (lon>275)&(lon<305);
%% SLHF
figure(1)
set(gcf,'position',[ 137          69        1131         735],'color','w')

subplot(2,2,year)
slhf(lsm>0) = NaN;
% imagesc(lon(GS_lon),lat(GS_lat),sshf_current(GS_lon,GS_lat)')
contourf(X(GS_lat,GS_lon)./1000,Y(GS_lat,GS_lon)./1000,slhf(GS_lon,GS_lat)',10)
title(['Surface Latent Heat Flux (DJFM)  ' datestr(time(1),'yyyy')])
set(gca,'ydir','reverse','xdir','reverse','fontsize',24)
xlabel('km')
ylabel('km')
axis equal
colorbar
subplot(2,2,[3 4])
plot(Y(GS_lat,find(GS_lon,1,'first'))./1000,nanmean(slhf(GS_lon,GS_lat)),'displayname',datestr(time(1),'yyyy'))
xlabel('km')
title_str = sprintf('Zonal Average [W m^{-2}]');
title(title_str)
set(gca,'fontsize',24)
hold on

%% SSHF
figure(2)
set(gcf,'position',[137   69  1131  735],'color','w')

subplot(2,2,year)
sshf(lsm>0) = NaN;
% imagesc(lon(GS_lon),lat(GS_lat),sshf_current(GS_lon,GS_lat)')
contourf(X(GS_lat,GS_lon)./1000,Y(GS_lat,GS_lon)./1000,sshf(GS_lon,GS_lat)',10)
title(['Surface Sensible Heat Flux (DJFM)  ' datestr(time(1),'yyyy')])
set(gca,'ydir','reverse','xdir','reverse','fontsize',24)
xlabel('km')
ylabel('km')
axis equal
colorbar
subplot(2,2,[3 4])
plot(Y(GS_lat,find(GS_lon,1,'first'))./1000,nanmean(sshf(GS_lon,GS_lat)),'displayname',datestr(time(1),'yyyy'))
xlabel('km')
title_str = sprintf('Zonal Average [W m^{-2}]');
title(title_str)
set(gca,'fontsize',24)
hold on


end

% contourf(X(GS_lat,GS_lon)./1000,Y(GS_lat,GS_lon)./1000,slhf(GS_lon,GS_lat)',10)
