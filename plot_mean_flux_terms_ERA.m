clear;clc;close all

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/lanczosfilter/')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

L = 500000; % m

intvl = 1; % look at every intvl'th timpepoint

plt_error_box = false;

year_vec = [2003:2019];

add_ssh_flag = false;
%
month_str = 'DJFM';

alpha_pos_flag = false;

filter_type = 'lanczos'; % filter type 'lanczos' or 'boxcar'

box_limits = [34 41; 143 168];

temporal_filter_length = 240;% hours (equiv to 10 days), min is 6



%% begin

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

temporal_mean_vec_L = zeros(length(year_vec),1);
temporal_mean_vec_S = zeros(length(year_vec),1);

cmap = cool(length(year_vec));

[n,m] = get_num_subplots( length(year_vec));
for i = 1:length(year_vec)
    year = year_vec(i);
    
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
        'time');
    
    load(sprintf('flux_terms_mean_ERA_%d_%d_%d_%d_%d',...
        box_limits(1),box_limits(3),box_limits(2),box_limits(4), year),...
        'spatial_mean_L','spatial_mean_S','temporal_mean_L','temporal_mean_S',...
        'temporal_filter_L','temporal_filter_S')
    
    % plot one number for each year
    figure(1)
    
    subplot(n,m,i)
    plot(time, temporal_filter_L,'-','linewidth',3,'displayname',sprintf('Latent %d',year))
    hold on
    plot(time, temporal_filter_S,'--','linewidth',3,'displayname',sprintf('Sensible %d',year))
    title(sprintf('filter size %d hours',temporal_filter_length),'interpreter','latex')
%     set(gca,'fontsize',20)
    xtickangle(-45)
%     legend('-dynamiclegend','location','best')
    
    set(gcf,'color','w','position',[ 232  1  1188  801],...
        'NumberTitle','off','Name',sprintf('flux terms'))
    temporal_mean_vec_L(i) = temporal_mean_L;
    temporal_mean_vec_S(i) = temporal_mean_S;
    
   
end
 update_figure_paper_size()
print(sprintf('imgs/flux_mean_ERA_time_filtered_%d_%d_%d_%d_L_%d_%s_%d_%d',...
    box_limits(1),box_limits(2),box_limits(3),box_limits(4),...
    L/1000,filter_type,year_vec(1),year_vec(end)),'-dpdf')

%%
load(sprintf('term68_%d_%s_%d',L/1000,filter_type,year_vec(1)),...
    'To_prime')
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon'); % other options 'lat','lon','patch_lat','patch_lon',
lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

figure(2)
subplot(1,2,1)
contourf(lon_er,lat_er,To_prime(:,:,1)')
hold on
y1 = box_limits(1);
y2 = box_limits(3);
x1 = box_limits(2);
x2 = box_limits(4);
plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'r','linewidth',3)

subplot(1,2,2)
year = year_vec(i);
% plot(year_vec,temporal_mean_vec_L./max(abs(temporal_mean_vec_L)),'displayname',sprintf('latent'),'linewidth',3)
plot(year_vec,temporal_mean_vec_L,'displayname',sprintf('latent'),'linewidth',3)
hold on
% plot(year_vec,temporal_mean_vec_S./max(abs(temporal_mean_vec_S)),'displayname',sprintf('sensible'),'linewidth',3)
plot(year_vec,temporal_mean_vec_S,'--','displayname',sprintf('sensible'),'linewidth',3)
legend('-dynamiclegend')
set(gca,'fontsize',20)
set(gcf,'color','w','position',[ 232  1  1188  801],...
    'NumberTitle','off','Name',sprintf('%s',filter_type))

update_figure_paper_size()
print(sprintf('imgs/flux_mean_ERA_%d_%d_%d_%d_L_%d_%s_%d_%d',...
    box_limits(1),box_limits(2),box_limits(3),box_limits(4),...
    L/1000,filter_type,year_vec(1), year_vec(end)),'-dpdf')

% show red box for T_o' or something?


