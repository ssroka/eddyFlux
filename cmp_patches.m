% good: SST d2m slhf sshf v10 u10 P0 t2m DT U_mag qa
% bad: qo qa

str = 'qo_patch';

DT_patch_old = load('/Volumes/SydneySroka_Remy/eddyFlux_data/global_ERA_data/ERA5_patch_data_2003.mat',str);
DT_patch_new = load('ERA5_patch_data_2003.mat',str);
subplot(1,3,1)
contourf(nanmean(eval(sprintf('DT_patch_old.%s',str)),3))
title('old')
colorbar
subplot(1,3,2)
contourf(nanmean(eval(sprintf('DT_patch_new.%s',str)),3))
title('new')
colorbar
subplot(1,3,3)
contourf(nanmean(eval(sprintf('DT_patch_new.%s',str)),3)-nanmean(eval(sprintf('DT_patch_old.%s',str)),3))
colorbar
title(str)
set(gcf,'position',[81         345        1255         460])