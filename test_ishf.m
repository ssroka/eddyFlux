% test instantaneous heat flux

clear;close all;clc
ishf_src = '/Volumes/SydneySroka_Anton/adaptor.mars.internal-1565102239.7902257-23416-1-12d8b204-b300-49a2-b629-002873207f66.nc';

ishf =  ncread(ishf_src,'ishf'); %
contourf(X(GS_lat,GS_lon)./1000,Y(GS_lat,GS_lon)./1000,sshf_GS',10)
colorbar
title(['SST [K] (DJFM)  ' datestr(time(1),'yyyy')])
set(gca,'ydir','reverse','xdir','reverse','fontsize',24)
ishf(lsm>0) = NaN;


imagesc(ishf')
figure

sshf = ncread(srcJFM,'sshf'); % J m^-2 / (24 hours in seconds) = W m^-2
S = sshf(:,:,1);
S(lsm>0) = NaN;

imagesc(S')


