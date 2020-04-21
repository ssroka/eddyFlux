ccc

filter_flag   = 'box'; % 'box' or 'zonal'

patch_str = 'Kur'; % 'GS'   'Kur'

% Kurishio
lat_bnds = [25 45];
lon_bnds = [130 170];
lat_sm_bnds = [30 40];

for year = 2003:2007

filename_data = sprintf('data_%s_%d_lat_%d_%d_lon_%d_%d.mat',filter_flag,year,lat_bnds(1),lat_bnds(2),lon_bnds(1),lon_bnds(2));
tic
load(filename_data,...
    'Lv','lat','lon','lat_bnds','lon_bnds','patch_lat','patch_lon',... % lat and lon vars
    'DT_patch','P0_patch','SST_patch','lsm_patch','slhf_patch','sshf_patch',... % ERA5 fields
    'u10_patch','v10_patch','time','t2m_patch','d2m_patch',...
    'U_mag','RH_patch','qa_patch','qo_patch')
load_time = toc;

tic
save(sprintf('/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/ERA5_patch_data_%d',year),...
    'Lv','lat','lon','lat_bnds','lon_bnds','patch_lat','patch_lon',... % lat and lon vars
    'DT_patch','P0_patch','SST_patch','lsm_patch','slhf_patch','sshf_patch',... % ERA5 fields
    'u10_patch','v10_patch','time','t2m_patch','d2m_patch',...
    'U_mag','RH_patch','qa_patch','qo_patch')           % calculated fields
save_time = toc;

fprintf('saved year %d, load time: %4.2f, save time: %4.2f \n',year,load_time,save_time)
    


end