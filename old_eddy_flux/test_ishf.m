% test instantaneous heat flux
% parameter conventions
% https://apps.ecmwf.int/codes/grib/param-db?id=146
% "The ECMWF convention for vertical fluxes is positive downwards."

clear;close all;clc


ishf_src = '/Volumes/SydneySroka_Anton/adaptor.mars.internal-1565102239.7902257-23416-1-12d8b204-b300-49a2-b629-002873207f66.nc';
ishf =  ncread(ishf_src,'ishf'); %

land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';
lsm = ncread(land_src,'lsm');

ishf(lsm>0) = NaN;
subplot(1,2,1)
imagesc(ishf')
colorbar

sshf =  ncread(ishf_src,'sshf'); %
sshf(lsm>0) = NaN;
subplot(1,2,2)
imagesc(sshf'/(60*60))
colorbar


% srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
% sshf = ncread(srcJFM,'sshf'); % J m^-2 / (24 hours in seconds) = W m^-2
% S = sshf(:,:,1);
% 
% S(lsm>0) = NaN;
% imagesc(S')

% % contourf(X(GS_lat,GS_lon)./1000,Y(GS_lat,GS_lon)./1000,sshf_GS',10)
% contourf(sshf_GS')
% colorbar
% title(['SST [K] (DJFM)  ' datestr(time(1),'yyyy')])
% set(gca,'ydir','reverse','xdir','reverse','fontsize',24)
% 
% figure




