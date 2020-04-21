
clear;close all;clc
% NOTE!: there was a sign change bug for the heat fluxes in case some plots
% made before 9/22 are not reproducable within a negative sign
addpath('~/Documents/MATLAB/util/')

flag_load_new_vars = false;
filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'box'; % 'box' or 'point'
err_box_lat = [32 38];
err_box_lon = [130 160];
error_units_flag = 'flux'; % 'flux' or 'power'

patch_str = 'Kur'; % 'GS'   'Kur'
land_src = '/Volumes/SydneySroka_Anton/ERA5_2018_Dec_1_0_LandSeaMask.nc';

alpha_eq_flag = 'constant_imposed'; % 'constant' (calculated from mean of flux (LS with constant valued model)
% 'linear'   (calculated from mean of flux (LS with model linear in T')
% 'constant_imposed' (pick the same constant for a_s and a_L
% sea level data from https://podaac.jpl.nasa.gov/dataset/RECON_SEA_LEVEL_OST_L4_V1
% ---------- Parameters --------------
% filter parameters
L = 500000; % [m] filtering box edge length

SSH_data = '/Volumes/SydneySroka_Remy/eddyFlux_data/CCAR_recon_sea_level_20000103_20090627_v1.nc';
time = ncread(SSH_data,'time')*24*3600;
time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');


for year = [3 7]
    
    switch year
        case 3 % 2003
            srcJFM = '/Volumes/SydneySroka_Anton/ERA5_2003_JFM_31d_6h.nc';
            srcD = '/Volumes/SydneySroka_Anton/ERA5_2002_Dec_31d_6h.nc';
            time_IDS = (time>datetime(2002,12,1)) & (time<datetime(2003,3,30));
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
            time_IDS = (time>datetime(2006,12,1)) & (time<datetime(2007,3,30));
            
    end
    
    %     time = double([ncread(srcD,'time');ncread(srcJFM,'time')])*60*60; % hours
    %     time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
    
    lat = ncread(SSH_data,'lat');
    lon = ncread(SSH_data,'lon');
    
    
    switch patch_str
        case 'GS'
            % Gulf Stream
            lat_bnds = [25 45];
            lon_bnds = [275 305];
        case 'Kur'
            % Kurishio
            lat_bnds = [25 45];
            lon_bnds = [130 170];
    end
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    contour(lon(patch_lon),lat(patch_lat),SSH(patch_lon,(patch_lat),i)',[10 10],'k')
    contourf(lon(patch_lon),lat(patch_lat),SSH_a(patch_lon,patch_lat,i)',[10],'k')
    
    
end

















