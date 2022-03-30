


for y = 2003:2018
   load(sprintf('ERA5_patch_data_%d.mat',y))
   save(sprintf('ERA5_SST_%d.mat',y),'SST_patch','lat','lon','patch_lat','patch_lon','lat_bnds','lon_bnds','lsm_patch')
end