% data from
%{
https://resources.marine.copernicus.eu/?option=com_csw&task=results
https://resources.marine.copernicus.eu/?option=com_csw&view=order&record_id=c0635fc4-07d3-4309-9d55-cfd3e6aa788b

%}
function [] = plot_SSH_contour(data_src,year,box_opt)
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
    'lat','lon','patch_lat','patch_lon');

% SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/CCAR_SSH_data/SSH_20000103_20090627_CCAR.mat','lon','lat','time','SSH');
    if year<2008 & year>2002
        SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/ssh_data_global/global-reanalysis-phy-001-030-daily_20021201_20070401.mat','lon','lat','time','SSH');
    elseif year<2013 & year>2007
        SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/ssh_data_global/global-reanalysis-phy-001-030-daily_20071201_20120331.mat','lon','lat','time','SSH');
    elseif year>2012 & year<=2018
        SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/ssh_data_global/global-reanalysis-phy-001-030-daily_20131201_20180331.mat','lon','lat','time','SSH');
    end
setup_lat_lon_vec;

patch_lat_bnds = [min(lat(patch_lat)) max(lat(patch_lat))];
patch_lon_bnds = [min(lon(patch_lon)) max(lon(patch_lon))];

patch_lat_SSH = (SSH_data.lat>=patch_lat_bnds(1)) & (SSH_data.lat<=patch_lat_bnds(2));
patch_lon_SSH = (SSH_data.lon>=patch_lon_bnds(1)) & (SSH_data.lon<=patch_lon_bnds(2));

        ssh_lvl = 0.5;

% switch year
%     case 2003
%         ssh_lvl = 0.5;
%     case 2007
%         ssh_lvl = 0.5;
% end

count = 1;
clear SSH_4_mean
for i = 1:length(SSH_data.time)
    Dec_SSH = (SSH_data.time.Year(i)==year-1) && (SSH_data.time.Month(i) == 12) && (mod(SSH_data.time.Day(i),7) == 0);
    JFM_SSH = (SSH_data.time.Year(i)==year) && (SSH_data.time.Month(i) <= 3) && (mod(SSH_data.time.Day(i),10) == 0);
    if (Dec_SSH || JFM_SSH) & (mod(i,2)==0)
        SSHlon=SSH_data.lon(patch_lon_SSH);
        SSHlat=SSH_data.lat(patch_lat_SSH);
        SSH_lat_box_TF = SSHlat>=box_opt(1,1) & SSHlat<=box_opt(1,2);
        SSH_lon_box_TF = SSHlon>=box_opt(2,1) & SSHlon<=box_opt(2,2);
        SSH_patch = SSH_data.SSH(patch_lon_SSH,patch_lat_SSH,i);
        [~,ch] = contour(SSHlon(SSH_lon_box_TF),SSHlat(SSH_lat_box_TF),...
            SSH_patch(SSH_lon_box_TF,SSH_lat_box_TF)',[1 1]*ssh_lvl);
        hold on
        set(ch,'linewidth',1,'color','k')
        count = count + 1;
    end
end


end

