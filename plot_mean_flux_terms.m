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


% nonlinear derivation terms
term_str_S = {...
    '\rho_a c_p C_D^s \Delta T U',...
    '\rho_a c_p C_D^s \Delta T \alpha_s T_o'' U''',...
    '\rho_a c_p C_D^s \Delta T'' \alpha_s T_o U''',...
    '\rho_a c_p C_D^s \Delta T'' U''',...
    '\rho_a c_p C_D^s \Delta T'' \alpha_s T_o''U''',...
    '\rho_a c_p C_D^s \Delta T \alpha_s T_o''U',...
    '\rho_a c_p C_D^s \Delta T \alpha_s T_o U''',...
    '\rho_a c_p C_D^s \Delta T'' \alpha_s T_o U',...
    };
%% begin

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

% create legend
leg_ent = cell(8,1);
for j = 1:8
    leg_ent = sprintf('term%d',j);
end

temporal_mean_vec = zeros(length(year_vec),8);

for i = 1:length(year_vec)
    year = year_vec(i);
    
    load(sprintf('flux_terms_mean_%d_%sfilt_%s_%d_%d_%d_%d_%d',L/1000,con_str,filter_type,...
        box_limits(1),box_limits(2),box_limits(3),box_limits(4), year),...
        'temporal_mean','temporal_filter') % other options 'spatial_mean',
    
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
        'time'); % other options 'lat','lon','patch_lat','patch_lon',
    
    % plot one number for each year
    figure(1)
    clf
    for j = 1:2
        if j == 1
            title_str = 'sensible';
        else
            title_str = 'latent';
        end
        
        for k = 1:8
            subplot(2,4,k)
            hold on
            plot(time, temporal_filter(:,k,j),'-','linewidth',3,'displayname',title_str)
            title(sprintf('$$%s [W m^{-2}]$$',term_str_S{k}),'interpreter','latex')
            set(gca,'fontsize',20)
            xtickangle(-45)
            legend('-dynamiclegend','location','best')
            ax = gca;
            ax.YAxis.Exponent = 0;
        end
        
    end
    set(gcf,'color','w','position',[43           1        1377         801],...
        'NumberTitle','off','Name',sprintf('flux terms %d',year))
    temporal_mean_vec(i,:) = sum(temporal_mean,3);
    
    update_figure_paper_size()
    print(sprintf('imgs/flux_terms_mean_sens_lat_%d_%d_%d_%d_L_%d_%s_%d',...
        box_limits(1),box_limits(2),box_limits(3),box_limits(4),...
        L/1000,filter_type,year),'-dpdf')
end

load(sprintf('term68_%d_%s_%d',L/1000,filter_type,year_vec(1)),...
    'To_prime')
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon'); % other options 'lat','lon','patch_lat','patch_lon',
lat_er = lat(patch_lat);
lon_er = lon(patch_lon);
%%
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
for i = 1:8
    plot(year_vec,temporal_mean_vec(:,i)./max(abs(temporal_mean_vec(:,i))),'displayname',sprintf('term %d',i),'linewidth',3)
    hold on
end
legend('-dynamiclegend')
set(gca,'fontsize',20)
set(gcf,'color','w','position',[ 232  1  1188  801],...
    'NumberTitle','off','Name',sprintf('%s %d',title_str,year))

update_figure_paper_size()
print(sprintf('imgs/flux_terms_mean_%d_%d_%d_%d_L_%d_%s_%d_%d',...
    box_limits(1),box_limits(2),box_limits(3),box_limits(4),...
    L/1000,filter_type,year_vec(1), year_vec(end)),'-dpdf')

% show red box for T_o' or something?


