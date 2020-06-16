clear;clc;close all

year_vec = [2003:2019];

data_src = '/Volumes/SydneySroka_Remy/eddyFlux_data/';

save_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

land_src = '/Volumes/SydneySroka_Remy/eddyFlux_data/LandSeaMask_20020101_00.nc';

% Kurishio
lat_bnds = [25 45];
lon_bnds = [130 170];

%% begin

lsm = ncread(land_src,'lsm');

for i = 1:length(year_vec)
    year = year_vec(i);
    
    srcD = sprintf('/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_%d_DJFM_31d_6h.nc',year-1);
    srcJFM = sprintf('/Volumes/SydneySroka_Remy/eddyFlux_data/ERA5_%d_DJFM_31d_6h.nc',year);
    
    time_DEC = double(ncread(srcD,'time'))*60*60; % hours
    time_DEC = datetime(time_DEC,'ConvertFrom','epochtime','epoch','1900-01-01');
    DEC_ids = time_DEC.Month==12;
    time_DEC = time_DEC(DEC_ids);
    
    time_JFM = double(ncread(srcJFM,'time'))*60*60; % hours
    time_JFM = datetime(time_JFM,'ConvertFrom','epochtime','epoch','1900-01-01');
    JFM_ids = time_JFM.Month==1 | time_JFM.Month==2 | time_JFM.Month==3;
    time_JFM = time_JFM(JFM_ids);
    
    time = [time_DEC;time_JFM];
    
    nt = length(time);
    
    lat = ncread(srcJFM,'latitude');
    lon = ncread(srcJFM,'longitude');
    
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    lsm_patch = lsm(patch_lon,patch_lat);
    
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [slhf_patch] = get_patch('slhf',srcD,DEC_ids,srcJFM,JFM_ids,...
        1./(-1*60*60),patch_lon,patch_lat,lsm_patch,true);
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [sshf_patch] = get_patch('sshf',srcD,DEC_ids,srcJFM,JFM_ids,...
        1./(-1*60*60),patch_lon,patch_lat,lsm_patch,true);
    % SST [K]
    [SST_patch] = get_patch('sst',srcD,DEC_ids,srcJFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % t2m [K]
    [t2m_patch] = get_patch('t2m',srcD,DEC_ids,srcJFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % d2m [K]
    [d2m_patch] = get_patch('d2m',srcD,DEC_ids,srcJFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % sp [Pa]
    [P0_patch] = get_patch('sp',srcD,DEC_ids,srcJFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % u10 [m/s]
    [u10_patch] = get_patch('u10',srcD,DEC_ids,srcJFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % v10 [m/s]
    [v10_patch] = get_patch('v10',srcD,DEC_ids,srcJFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
        
    
    DT_patch = SST_patch - t2m_patch;
    U_mag = sqrt(u10_patch.^2+v10_patch.^2);
    TDC = d2m_patch-273.15;
    TC  = t2m_patch-273.15;
    e0 = 0.6112*exp((17.67.*TDC)./(243.5+TDC))*1000; % Pa
    e_sat = SAM_psatWater(t2m_patch); % Pa
    RH_patch = 100.*e0./e_sat;
    
    % ---- qo ----
    qo_patch = SAM_qsatWater(SST_patch, P0_patch) ;
    lsm_patch_mat = repmat(lsm_patch,1,1,nt);
    qo_patch(lsm_patch_mat>0)=NaN;
    
    qa_patch = zeros(size(v10_patch,1),size(v10_patch,2),nt);
    for tt=1:nt
        % ---- qa ----
%         e_sat = SAM_psatWater(t2m_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat(:,:,tt);
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch(:,:,tt) = r./(1+r);
    end
    
 
    
    save(sprintf('ERA5_patch_data_%d',year),'slhf_patch','sshf_patch','SST_patch',...
        't2m_patch','d2m_patch','P0_patch','P0_patch','u10_patch','v10_patch',...
        'qo_patch','qa_patch','DT_patch','RH_patch','U_mag',...
        'lat','lat_bnds','lon','lon_bnds','lsm_patch','patch_lat','patch_lon',...
        'time')
    
    
end

