
data_base = '/Volumes/SydneySroka_Remy/eddyFlux/';
data_src = [data_base 'ERA5_data/'];

year_vec = [2003:2018];
box_num = 9;
er_box_num = 3;

switch box_num
    case 1
        box_opt = [36 41.5; 143 152];
    case 2
        box_opt = [30 44.5; 148 169];
    case 3
        box_opt = [30 41.5; 142.5 169];
    case 4
        box_opt = [30 41.5; 142.5 153];
    case 5
        box_opt = [30 41.5; 153 169];
    case 6
        box_opt = [34 38; 142.5 169];
    case 7
        box_opt = [30 33.83; 142.5 169];
    case 8
        box_opt = [33.83 37.663; 142.5 169];
    case 9
        box_opt = [37.663 41.5; 142.5 169];
end

switch er_box_num
    case 3
        er_box = [35.75 41.5; 142.5 155.7500];% [30 41.5; 142.5 169];%
end

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');
setup_lat_lon_vec;
HF = zeros(1,length(year_vec));
for j = 1:length(year_vec)
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(j)),...
        'sshf_patch','slhf_patch');
    
    SSHF_mean = nanmean(nanmean(nanmean(sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),3)));
    SLHF_mean = nanmean(nanmean(nanmean(slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),3)));
    
    HF(j) =   SSHF_mean+  SLHF_mean;
    
end

save(sprintf('HF_3_boxes_%d',box_num),'HF')















