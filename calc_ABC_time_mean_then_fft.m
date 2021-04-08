clear;clc;close all

year_vec = [2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/lanczosfilter/')
addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/fft_filter/')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

L = 250000; % m

intvl = 1; % look at every intvl'th timpepoint

alpha_pos_flag = false;

fft_first_flag = false;

filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'

%% filter set up

load('env_const.mat')

files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),'SST_patch','lat','lon','patch_lat','patch_lon');


d_lat = abs(files_for_size.lat(2)-files_for_size.lat(1));
d_lon = abs(files_for_size.lon(2)-files_for_size.lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
Ny = floor(L/dy)+mod(floor(L/dx),2)+1; % make Ny odd

NaN_inds = isnan(files_for_size.SST_patch(:,:,1));

if strcmp(filter_type,'fft')
    cf = (1/(2*L));
    box_opt = [32 38; 143 167];
elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
end

% this is the same for every year
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');

file_for_box = load(sprintf('Qs_QL_optimization_data_L_%d_filt_%s_%d',L/1000,filter_type,2003),'box_limits');

box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

% to index out of *_patch fields
opt_patch_lat = box_lat(patch_lat);
opt_patch_lon = box_lon(patch_lon);

prime_lat =  lat>=file_for_box.box_limits(1,1) & lat<=file_for_box.box_limits(1,2);
prime_lon = lon>=file_for_box.box_limits(2,1) & lon<=file_for_box.box_limits(2,2);

% to index out of *_prime fields
opt_prime_lat = box_lat(prime_lat);
opt_prime_lon = box_lon(prime_lon);

lat_plot = lat(box_lat);
lon_plot = lon(box_lon);

m = sum(opt_prime_lon);
n = sum(opt_prime_lat);

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

salinity = 34*ones(m,n);% ppt for Lv calculation

if strcmp(filter_type,'boxcar')
    [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
end

if fft_first_flag
    fft_str = '';
else
    fft_str = 'noFFT_';
end


for i = 1:length(year_vec)
    year = year_vec(i);
    
    load(sprintf('ABC_terms_%d_%sfilt_%s%s_%d',L/1000,con_str,fft_str,filter_type,year))
    
    figure(1)
    ax = subplot(2,4,1);
    A_time_mean = FFT2D_filter(nanmean(A,3),dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,A_time_mean');
    title('A:$\rho_a C_D\overline{U}\overline{\Delta h}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(B,3)');
    title('B:$\rho_a C_D\overline{U''\Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C1,3)');
    title('C1:$\rho_a C_D\overline{T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,4);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C2,3)');
    title('C2:$\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,5);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C3,3)');
    title('C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,6);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(D,3)');
    title('D:$\rho_a C_D\overline{U''\overline{\Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,7);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(E1,3)');
    title('E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,8);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(E2,3)');
    title('E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    set(gcf,'color','w','position',[61 221 1359 581],'NumberTitle','off','Name',num2str(year))
    update_figure_paper_size()
    print(sprintf('imgs/ABC_L_%d_%s%s_%d',L/1000,fft_str,filter_type,year),'-dpdf')
    %%
    figure(2)
    
    ax = subplot(2,2,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(A+B+C1+C2+C3+D+E1+E2,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    set(gca,'clim',[200 400])
    
    ax = subplot(2,2,2);
    [~,h] = contourf(lon(opt_patch_lon), lat(opt_patch_lat),...
        nanmean(slhf_patch(opt_patch_lon,opt_patch_lat,:)+sshf_patch(opt_patch_lon,opt_patch_lat,:),3)');
    title('SLHF + SSHF from ERA5 Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    set(gca,'clim',[200 400])
    
    ax = subplot(2,2,3);
    [~,h] = contourf(lon(opt_patch_lon), lat(opt_patch_lat),...
        nanmean(slhf_patch_bar+sshf_patch_bar,3)');
    title('$\overline{SSHF + SLHF}$ from ERA5 Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    set(gca,'clim',[200 400])
    
    figure(3)
    ax = subplot(4,2,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(dh_CTRL,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(4,2,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(dh_prime,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(4,2,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(adh_CTRL,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(4,2,4);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(adh_prime,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    
    ax = subplot(4,2,5);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(U_CTRL,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(4,2,6);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(U_prime,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    
    ax = subplot(4,2,7);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(To_prime,3)');
    title('total:$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
end


function [] = format_fig(h,plt_num,max_val,min_val)

set(h,'edgecolor','none')
set(gca,'ydir','normal','fontsize',15)
colorbar
xlabel('deg')
ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end


% get_eddy_contribution_plot;
