clear;close all;clc
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
patch_str = 'Kur'; % 'GS'   'Kur'
% land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';


dLat = 0.25;
Lat2m = 111000;
dLon = 0.25;
Lon2m = 85400;

% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length
clrs = 'kkkkkkkkkk';
for yy = [3:7]
    
    switch yy
        case 3 % 2003
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
        case 4 % 2004
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2004_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2003_Dec_31d_6h.nc';
        case 5 % 2005
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2005_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2004_Dec_31d_6h.nc';
        case 6 % 2006
            srcJFM = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2006_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_2005_Dec_31d_6h.nc';
        case 7 % 2007
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2007_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2006_Dec_31d_6h.nc';
    end
    
    %     time = double([ncread(srcD,'time');ncread(srcJFM,'time')])*60*60; % hours
    %     time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
    
    %     lat = ncread(srcJFM,'latitude');
    %     lon = ncread(srcJFM,'longitude');
    time.Year = yy+2000;
    switch patch_str
        case 'GS'
            % Gulf Stream
            lat_bnds = [25 45];
            lon_bnds = [275 305];
        case 'Kur'
            % Kurishio
            lat_bnds = [25 45];
            lon_bnds = [130 170];
            lat_sm_bnds = [30 40];
    end
    
    filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,time(end).Year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
    load(filename_data)
    
    g = var(SST_prime,0,3);
    V(yy) =nanmedian(g(:));
    patch_lat_sm = (lat>lat_sm_bnds(1))&(lat<lat_sm_bnds(2));
    inds_sm = find(patch_lat_sm(patch_lat));
    HF = sshf_patch(:,inds_sm,:)+slhf_patch(:,inds_sm,:); % Watts m^-2
    H(yy) = nansum(nansum(nansum(HF*dLat*Lat2m*dLon*Lon2m))); % Watts
    
    figure(1)
    subplot(2,5,(yy-2)*(yy<6)+(yy>5)*yy)
    contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(HF,3)',30,'k')
    colorbar
    set(gca,'clim',[150 550]);
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    title(['mean flux in ' num2str(time(end).Year)])
    
    
    avgHF(yy) = nanmedian(HF(:));
    T(yy) = time(end).Year;
    s_str = sprintf('SSHF [W/$$m^{-2}$$] (mean DJFM %d)',time(end).Year);
    L_str = sprintf('SLHF [W/$$m^{-2}$$] (mean DJFM %d)',time(end).Year);
    
    subplot(2,5,[4,5,9,10])
%     plot(V(yy),H(yy)/H(3),'*','displayname',num2str(time(end).Year))
    plot(V(yy),avgHF(yy)/avgHF(3),'*','displayname',num2str(time(end).Year))
    hold on
    
end
title(sprintf('ERA5 SSHF  (median DJFM %d)',time(end).Year),'interpreter','latex')
xlabel('variance')
lh = legend('-dynamiclegend');
set(lh,'interpreter','latex','location','best')
set(gca,'fontsize',20,'YAxisLocation','right')
set(gcf,'color','w','position',[ 56         264        1285         541])

% regression
Y = avgHF(3:7)'/avgHF(3);
% Y = H(3:7)'./H(3);
X = [ones(length(3:7),1) V(3:7)'];
[B,~,~,~,STAT] = regress(Y,X);
p = polyfit(X(:,2),Y,1);
v_vec = [min(V(3:7)) max(V(3:7))];
% ylabel('normalized by total integrated flux in 2003')
ylabel('normalized by flux in 2003')
plot(v_vec,polyval(p,v_vec),'r-','linewidth',2,'displayname','linear fit')
% title(['$$\\V$$ as a function of total flux $$[W]$$, $$R^2=$$' num2str(STAT(1))] ,'interpreter','latex')
title(['$$\\V$$ as a function of total flux $$[W m^{-2}]$$, $$R^2=$$' num2str(STAT(1))] ,'interpreter','latex')
% addpath('~/Documents/MATLAB/util/tightfig/') %[W/$$m^{-2}$$]
% tightfig(gcf)
set(gcf,'color','w','position',[ 56         264        1285         541])



%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/var_vs_flux','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/var_vs_flux_SST_prime','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/var_vs_flux_SST_prime_med','-dpng')


%}


