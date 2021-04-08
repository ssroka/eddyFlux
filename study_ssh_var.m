% data from
%{
https://resources.marine.copernicus.eu/?option=com_csw&task=results
https://resources.marine.copernicus.eu/?option=com_csw&view=order&record_id=c0635fc4-07d3-4309-9d55-cfd3e6aa788b

%}
clear;clc;close all

year_vec = [2003];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

%     box_opt = [32 38; 143 167];
%     box_opt = [30 42; 144 168];
box_opt = [30 44.5; 148 169];

%% filter set up
load('env_const.mat')

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'SST_patch','lat','lon','patch_lat','patch_lon',...
    'slhf_patch','sshf_patch');

var_SSH = zeros(length(year_vec),1);

for j = 1:length(year_vec)
    year = year_vec(j);
    
    fprintf('year %d\n',year)
    
    % SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/CCAR_SSH_data/SSH_20000103_20090627_CCAR.mat','lon','lat','time','SSH');
    if year<2008 & year>2002
        SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/ssh_data_global/global-reanalysis-phy-001-030-daily_20021201_20070401.mat','lon','lat','time','SSH');
    elseif year<2013 & year>2007
        SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/ssh_data_global/global-reanalysis-phy-001-030-daily_20071201_20120331.mat','lon','lat','time','SSH');
    elseif year>2012 & year<2018
        SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/ssh_data_global/global-reanalysis-phy-001-030-daily_20131201_20180331.mat','lon','lat','time','SSH');
        
    end
    
%     patch_lat_bnds = [min(lat(patch_lat)) max(lat(patch_lat))];
%     patch_lon_bnds = [min(lon(patch_lon)) max(lon(patch_lon))];
    
    patch_lat_bnds = box_opt(1,1:2);
    patch_lon_bnds = box_opt(2,1:2);
    
    patch_lat_SSH = (SSH_data.lat>=patch_lat_bnds(1)) & (SSH_data.lat<=patch_lat_bnds(2));
    patch_lon_SSH = (SSH_data.lon>=patch_lon_bnds(1)) & (SSH_data.lon<=patch_lon_bnds(2));
        
    % switch year
    %     case 2003
    %         ssh_lvl = 0.5;
    %     case 2007
    %         ssh_lvl = 0.5;
    % end
    SSH = zeros(sum(patch_lat_SSH),sum(patch_lon_SSH),length(SSH_data.time));
    count = 1;
    clear SSH_4_mean
    for i = 1:length(SSH_data.time)
        Dec_SSH = (SSH_data.time.Year(i)==year-1) && (SSH_data.time.Month(i) == 12);
        JFM_SSH = (SSH_data.time.Year(i)==year) && (SSH_data.time.Month(i) <= 3);
        if (Dec_SSH || JFM_SSH)
            SSH(:,:,count) = SSH_data.SSH(patch_lon_SSH,patch_lat_SSH,i)';
            count = count + 1;
        end
    end
    SSH = SSH(:,:,1:count-1);
    figure(year)
    [~,ch] = contourf(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        var(SSH,0,3));
    clim = get(gca,'clim');
    hold on
    [C,h] = contour(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        mean(SSH,3),'w','linewidth',2);
    clabel(C,h,'color','w')
    set(gca,'clim',clim)
    colorbar
    var_SSH(j) = mean(mean(var(SSH,0,3)));
    set(gca,'fontsize',25)
    set(gcf,'color','w','position',[1 215 1439 587],'NumberTitle','off','Name',num2str(year))
    title(['SSH variance (colorbar), white contours = time mean',num2str(year)])
    update_figure_paper_size()
    print(sprintf('imgs/var_ssh_%d',year),'-dpdf')
 
end

save(sprintf('var_SSH_%d_%d',year_vec(1),year_vec(end)),'var_SSH','year_vec')



