clear;clc;close all

year_vec = [2003];


% data_src = '/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/';
% save_src = '/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/';
data_src = '~/Downloads/NCEP_data/';
save_src = '~/Downloads/NCEP_data/';

addpath('/Volumes/SydneySroka_Remy/eddyFlux/filter')
addpath('/Volumes/SydneySroka_Remy/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

% land_src = '/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/land.nc';
land_src = ['/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/lsm_NCEP.nc'];

% Kurishio
lat_bnds = [25 45];
lon_bnds = [130 170];

%% begin

lsm = double(ncread(land_src,'LAND'));

for i = 1:length(year_vec)
    year = year_vec(i);
    
    Ta = sprintf('%sT2m_NCEP_%d.nc',data_src,year);
    
    SST = sprintf('%sSST_NCEP_%d.nc',data_src,year);
    
    slhf = sprintf('%sslhf_NCEP_%d.nc',data_src,year);
    
    sshf = sprintf('%ssshf_NCEP_%d.nc',data_src,year);
    
    u10 = sprintf('%su10_NCEP_%d.nc',data_src,year);
    v10 = sprintf('%sv10_NCEP_%d.nc',data_src,year);
    
    P0 = sprintf('%sp0_NCEP_%d.nc',data_src,year);
    
    qa = sprintf('%sqa_NCEP_%d.nc',data_src,year); % [kg/kg]

    
    time = double(ncread(Ta,'T'))*24*60*60; % hours to seconds
    time = datetime(time,'ConvertFrom','epochtime','epoch','1948-01-01 12:00:00');
    
    nt = length(time);
    
    lat = double(ncread(Ta,'Y'));
    lon = double(ncread(Ta,'X'));
    
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    lsm_patch = lsm(patch_lon,patch_lat,1);
    
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [slhf_patch] = get_patch_NCEP('heat_flux',slhf,...
        patch_lon,patch_lat,lsm_patch,true);
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive d    qownward" convention
    [sshf_patch] = get_patch_NCEP('heat_flux',sshf,...
        patch_lon,patch_lat,lsm_patch,true);
    % SST [K]
    [SST_patch] = get_patch_NCEP('temp',SST,...
        patch_lon,patch_lat,lsm_patch,true);
    % t2m [K]
    [t2m_patch] = get_patch_NCEP('temp',Ta,...
        patch_lon,patch_lat,lsm_patch,true);
    % qa [%]
    [qa_patch] = get_patch_NCEP('qa',qa,...
        patch_lon,patch_lat,lsm_patch,true);
    % sp [Pa]
    [P0_patch] = get_patch_NCEP('pressure',P0,...
        patch_lon,patch_lat,lsm_patch,true);
    % u10 [m/s]
    [u10_patch] = get_patch_NCEP('u',u10,...
        patch_lon,patch_lat,lsm_patch,true);
    % v10 [m/s]
    [v10_patch] = get_patch_NCEP('v',v10,...
        patch_lon,patch_lat,lsm_patch,true);
        
    
    DT_patch = SST_patch - t2m_patch;
    U_mag = sqrt(u10_patch.^2+v10_patch.^2);
%     TDC = d2m_patch-273.15;
%     TC  = t2m_patch-273.15;
%     e0 = 0.6112*exp((17.67.*TDC)./(243.5+TDC))*1000; % Pa
    e_sat = SAM_psatWater(t2m_patch); % Pa
%     RH_patch = 100.*e0./e_sat;
    
    % ---- qo ----
    qo_patch = SAM_qsatWater(SST_patch, P0_patch) ;
    lsm_patch_mat = repmat(lsm_patch,1,1,nt);
    qo_patch(lsm_patch_mat>0)=NaN;
    
%     qa_patch = zeros(size(v10_patch,1),size(v10_patch,2),nt);
%     for tt=1:nt
%         % ---- qa ----
% %         e_sat = SAM_psatWater(t2m_patch(:,:,tt));
%         e = RH_patch(:,:,tt)./100.*e_sat(:,:,tt);
%         r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
%         qa_patch(:,:,tt) = r./(1+r);
%     end
    
 
    
    save(sprintf('NCEP_patch_data_%d',year),'slhf_patch','sshf_patch','SST_patch',...
        't2m_patch','P0_patch','u10_patch','v10_patch',...
        'qo_patch','qa_patch','DT_patch','U_mag',...
        'lat','lat_bnds','lon','lon_bnds','lsm_patch','patch_lat','patch_lon',...
        'time')
    
    
end

