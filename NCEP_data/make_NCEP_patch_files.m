clear;clc;close all

year_vec = [2003];


data_src = '/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/';
save_src = '/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/';

addpath('/Volumes/SydneySroka_Remy/eddyFlux/filter')
addpath('/Volumes/SydneySroka_Remy/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

land_src = '/Volumes/SydneySroka_Remy/eddyFlux/NCEP_data/land.nc';

% Kurishio
lat_bnds = [25 45];
lon_bnds = [130 170];

%% begin

lsm = ncread(land_src,'land');

for i = 1:length(year_vec)
    year = year_vec(i);
    
    Ta_D = sprintf('%sair_2m_gauss_%d.nc',data_src,year-1);
    Ta_JFM = sprintf('%sair_2m_gauss_%d.nc',data_src,year);
    
    SST_D = sprintf('%sskt_sfc_gauss_%d.nc',data_src,year-1);
    SST_JFM = sprintf('%sskt_sfc_gauss_%d.nc',data_src,year);
    
    slhf_D = sprintf('%slhtfl_sfc_gauss_%d.nc',data_src,year-1);
    slhf_JFM = sprintf('%slhtfl_sfc_gauss_%d.nc',data_src,year);
    
    sshf_D = sprintf('%sshtfl_sfc_gauss_%d.nc',data_src,year-1);
    sshf_JFM = sprintf('%sshtfl_sfc_gauss_%d.nc',data_src,year);
    
    u10_D = sprintf('%suwnd_10m_gauss_%d.nc',data_src,year-1);
    u10_JFM = sprintf('%suwnd_10m_gauss_%d.nc',data_src,year);
    
    v10_D = sprintf('%svwnd_10m_gauss_%d.nc',data_src,year-1);
    v10_JFM = sprintf('%svwnd_10m_gauss_%d.nc',data_src,year);
    
    P0_D = sprintf('%sslp_%d.nc',data_src,year-1);
    P0_JFM = sprintf('%sslp_%d.nc',data_src,year);
    
    RH_D = sprintf('%srhum_sig995_%d.nc',data_src,year-1);
    RH_JFM = sprintf('%srhum_sig995_%d.nc',data_src,year);

    
    
    time_DEC = double(ncread(Ta_D,'time'))*60*60; % hours to seconds
    time_DEC = datetime(time_DEC,'ConvertFrom','epochtime','epoch','1800-01-01');
    DEC_ids = time_DEC.Month==12;
    time_DEC = time_DEC(DEC_ids);
    
    time_JFM = double(ncread(Ta_JFM,'time'))*60*60; % hours
    time_JFM = datetime(time_JFM,'ConvertFrom','epochtime','epoch','1800-01-01');
    JFM_ids = time_JFM.Month==1 | time_JFM.Month==2 | time_JFM.Month==3;
    time_JFM = time_JFM(JFM_ids);
    
    time = [time_DEC;time_JFM];
    
    nt = length(time);
    
    lat = double(ncread(Ta_JFM,'lat'));
    lon = double(ncread(Ta_JFM,'lon'));
    
    patch_lat = (lat>lat_bnds(1))&(lat<lat_bnds(2));
    patch_lon = (lon>lon_bnds(1))&(lon<lon_bnds(2));
    
    lsm_patch = lsm(patch_lon,patch_lat);
    
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [slhf_patch] = get_patch_NCEP('lhtfl',slhf_D,DEC_ids,slhf_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % SSHF J m^-2 / (24 hours in seconds) = W m^-2 and (-1) since ERA5
    % has a "positive downward" convention
    [sshf_patch] = get_patch_NCEP('shtfl',sshf_D,DEC_ids,sshf_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % SST [K]
    [SST_patch] = get_patch_NCEP('skt',SST_D,DEC_ids,SST_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % t2m [K]
    [t2m_patch] = get_patch_NCEP('air',Ta_D,DEC_ids,Ta_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % RH [%]
    [RH_patch] = get_patch_NCEP('rhum',RH_D,DEC_ids,RH_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % sp [Pa]
    [P0_patch] = get_patch_NCEP('slp',P0_D,DEC_ids,P0_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % u10 [m/s]
    [u10_patch] = get_patch_NCEP('uwnd',u10_D,DEC_ids,u10_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
    % v10 [m/s]
    [v10_patch] = get_patch_NCEP('vwnd',v10_D,DEC_ids,v10_JFM,JFM_ids,...
        1,patch_lon,patch_lat,lsm_patch,true);
        
    
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
    
    qa_patch = zeros(size(v10_patch,1),size(v10_patch,2),nt);
    for tt=1:nt
        % ---- qa ----
%         e_sat = SAM_psatWater(t2m_patch(:,:,tt));
        e = RH_patch(:,:,tt)./100.*e_sat(:,:,tt);
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch(:,:,tt) = r./(1+r);
    end
    
 
    
    save(sprintf('NCEP_patch_data_%d',year),'slhf_patch','sshf_patch','SST_patch',...
        't2m_patch','P0_patch','u10_patch','v10_patch',...
        'qo_patch','qa_patch','DT_patch','RH_patch','U_mag',...
        'lat','lat_bnds','lon','lon_bnds','lsm_patch','patch_lat','patch_lon',...
        'time')
    
    
end

