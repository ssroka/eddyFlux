
% SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/CCAR_SSH_data/SSH_20000103_20090627_CCAR.mat','lon','lat','time','SSH');
SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/global-reanalysis-phy-001-030-daily.mat','lon','lat','time','SSH');


patch_lat_bnds = [min(lat(patch_lat)) max(lat(patch_lat))];
patch_lon_bnds = [min(lon(patch_lon)) max(lon(patch_lon))];

patch_lat_SSH = (SSH_data.lat>=patch_lat_bnds(1)) & (SSH_data.lat<=patch_lat_bnds(2));
patch_lon_SSH = (SSH_data.lon>=patch_lon_bnds(1)) & (SSH_data.lon<=patch_lon_bnds(2));

switch year
    case 2003
        ssh_lvl = 80;
    case 2007
        ssh_lvl = 10;
end

count = 1;
clear SSH_4_mean
for i = 1:length(SSH_data.time)
    Dec_SSH = (SSH_data.time.Year(i)==year-1) && (SSH_data.time.Month(i) == 12);
    JFM_SSH = (SSH_data.time.Year(i)==year) && (SSH_data.time.Month(i) <= 3);
    if Dec_SSH || JFM_SSH
        SSH_4_mean(:,:,count) = SSH_data.SSH(patch_lon_SSH,patch_lat_SSH,i);
        [~,ch] = contour(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
            SSH_data.SSH(patch_lon_SSH,patch_lat_SSH,i)',[1 1]*ssh_lvl);
        hold on
        set(ch,'linewidth',1,'color','k')
        count = count + 1;
    end
end

% [~,ch] = contour(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
%     nanmean(SSH_4_mean,3)');
% set(ch,'linewidth',1,'color','k')
