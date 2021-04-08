


c = [-1000 1000];


[Ta_CTRL,Ta_prime] = FFT2D_filter(t2m_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
[qa_CTRL,qa_prime] = FFT2D_filter(qa_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);


subplot(2,2,1)
t1 = c_p_air.*SST_prime.*(SST_patch_CTRL-Ta_CTRL);
[t1_CTRL] = FFT2D_filter(t1,dx,cf,debug_flag,lat_plot,lon_plot);
[~,h] = contourf(lon_plot,lat_plot,t1_CTRL','linestyle','none');
colorbar
set(gca,'clim',c)

subplot(2,2,2)
t2 = c_p_air.*SST_prime.*(SST_prime-Ta_prime);
[t2_CTRL] = FFT2D_filter(t2,dx,cf,debug_flag,lat_plot,lon_plot);
[~,h] = contourf(lon_plot,lat_plot,t2_CTRL','linestyle','none');
colorbar
set(gca,'clim',c)

subplot(2,2,3)
t3 = Lv.*SST_prime.*(qo_CTRL-qa_CTRL);
[t3_CTRL] = FFT2D_filter(t3,dx,cf,debug_flag,lat_plot,lon_plot);
[~,h] = contourf(lon_plot,lat_plot,t3_CTRL','linestyle','none');
colorbar
set(gca,'clim',c)

subplot(2,2,4)
t4 = c_p_air.*SST_prime.*(qo_prime-qa_prime);
[t4_CTRL] = FFT2D_filter(t4,dx,cf,debug_flag,lat_plot,lon_plot);
[~,h] = contourf(lon_plot,lat_plot,t4_CTRL','linestyle','none');
colorbar
set(gca,'clim',c)

figure
t4 = SST_prime.*(qo_prime-qa_prime);
[~,h] = contourf(lon_plot,lat_plot,SST_prime','linestyle','none');
colorbar


