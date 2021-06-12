% data from
%{
https://resources.marine.copernicus.eu/?option=com_csw&task=results
https://resources.marine.copernicus.eu/?option=com_csw&view=order&record_id=c0635fc4-07d3-4309-9d55-cfd3e6aa788b

%}
% clear;clc;close all
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');

% filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year);

file_for_box = load(sprintf('beta_Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year_vec(1)),'box_opt');

box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

% to index out of *_patch fields
opt_patch_lat = box_lat(patch_lat);
opt_patch_lon = box_lon(patch_lon);

prime_lat =  lat>=file_for_box.box_opt(1,1) & lat<=file_for_box.box_opt(1,2);
prime_lon = lon>=file_for_box.box_opt(2,1) & lon<=file_for_box.box_opt(2,2);

% to index out of *_prime fields
opt_prime_lat = box_lat(prime_lat);
opt_prime_lon = box_lon(prime_lon);

lat_plot = lat(box_lat);
lon_plot = lon(box_lon);

m = sum(opt_prime_lon);
n = sum(opt_prime_lat);

if beta_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

salinity = 34*ones(m,n);% ppt for Lv calculation

if strcmp(filter_type,'boxcar')
    [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
end


year_vec = [2003 2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

%     box_opt = [32 38; 143 167];
%     box_opt = [30 42; 144 168];
% box_opt = [30 44.5; 148 169];
box_opt = [30 41.5; 142.5 169];
cntr = [0:0.2:0.8];

%% filter set up
load('env_const.mat')

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'SST_patch','lat','lon','patch_lat','patch_lon',...
    'slhf_patch','sshf_patch');

cntr_length = zeros(length(cntr),length(year_vec));
d = zeros(length(cntr),length(year_vec));

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
    
    meanSSH = mean(SSH,3);
    
    [~,ch] = contourf(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        meanSSH);
    clim = get(gca,'clim');
    hold on
    [C,h] = contour(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        meanSSH,cntr,'w','linewidth',2);
    clabel(C,h,'color','w')
    set(h,'levellist',cntr)
    set(gca,'clim',clim)
    colorbar
    set(gca,'fontsize',25)
    set(gcf,'color','w','position',[1 215 1439 587],'NumberTitle','off','Name',num2str(year))
    title(['SSH variance (colorbar), white contours = time mean',num2str(year)])
    C(:,C(1,:)<20) = [];
    Vq = interp2(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        meanSSH,C(1,:),C(2,:));
    for i = 1:length(cntr)
%         bM = regionprops(meanSSH>=cntr(i),'Perimeter');
%         cntr_length(i,j) = bM.Perimeter;
        inds = find(abs(Vq-cntr(i))<1e-10);
        for k = 1:length(inds)-1
            delta = norm(C(:,k+1)-C(:,k));
            if delta < 0.25
                d(i,j) = d(i,j) + delta;
            end
        end
    end

        load(sprintf('AbBbCb_terms_%d_%sfilt_%s%s_box%d_%d',L/1000,con_str,fft_str,filter_type,box_num,year))

      SSHF_mean = nanmean(nanmean(sshf_patch(opt_patch_lon,opt_patch_lat,1)));
    SLHF_mean = nanmean(nanmean(slhf_patch(opt_patch_lon,opt_patch_lat,1)));
    
HF(j) =   SSHF_mean+  SLHF_mean;
    %{
        inds = C(1,:)<20;
        C(:,inds) = [];
        plot(C(1,:),C(2,:),'mo')
        
        
        %}
%     d = 0;
%     for i = 1121:1454
%         d = d + norm(C(:,i+1)-C(:,i));
%     end
%         d2 = 0;
%     for i = 1098:1119
%         d2 = d2 + norm(C(:,i+1)-C(:,i));
%     end
%         d3 = 0;
%     for i = 1456:1543
%         d3 = d3 + norm(C(:,i+1)-C(:,i));
%     end
    
    %     update_figure_paper_size()
    %     print(sprintf('imgs/ssh_%d',year),'-dpdf')
    
end
d
figure
for i = 1:size(d,2)
    plot(cntr, d(:,i),'o-','linewidth',2,'displayname',num2str(year_vec(i)))
    hold on
end
legend('-dynamiclegend')
xlabel('contour height')
ylabel('distance [deg]')
set(gca,'fontsize',25,'xtick',cntr)
set(gcf,'color','w','position',[1 215 1439 587])

figure

plot(d(3,:),HF)




update_figure_paper_size()
print(sprintf('imgs/ssh_length_%d',year),'-dpdf')
    