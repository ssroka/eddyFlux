clear;clc;close all

year_vec = [2003];

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

filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'


param_3_str = '_3param';

fft_first_flag = true; % Do not change, USE calc_ABC_time_mean_then_fft if you want this
% option to be false

box_num = 2;

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
    %     box_opt = [32 38; 143 167];
    % box_opt = [30 42; 144 168];
%     box_opt = [30 44.5; 148 169];
%         box_opt = [36 41.5; 143 152];

elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
end

switch box_num
    case 1
        box_opt = [36 41.5; 143 152];
    case 2
        box_opt = [30 44.5; 148 169];
end

% this is the same for every year
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');

% filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year);

file_for_box = load(sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,2003),'box_opt');

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
    % fft_str = 'noFFT_'; USE calc_ABC_time_mean_then_fft if you want this
    % option
end

if param_3_str
    param_num_str = '_3param';
else
    param_num_str = '';
    
end

for i = 1:length(year_vec)
    year = year_vec(i);
    
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile,'qo_patch','qa_patch','SST_patch','DT_patch')
    load(sprintf('ABC_terms_%d_%sfilt_%s%s%s_box%d_%d_synthetic',L/1000,con_str,fft_str,filter_type,param_num_str,box_num,year))
    
    figure(1)
    ax = subplot(2,4,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(A,3)');
    title('A:$\rho_a C_D\overline{U}\overline{\Delta h}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(B,3)');
    title('B:$\rho_a C_D\overline{U''\Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C1,3)');
    title('C1:$\rho_a C_D\overline{U}\overline{T_o''\alpha \Delta h''}$','interpreter','latex')
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
    
    for i = 2:8
        subplot(2,4,i)
        set(gca,'clim',[-1 1])
    end
    set(gcf,'color','w','position',[61 221 1359 581])
    update_figure_paper_size()
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_L_%d_%s_box%d_%d_synthetic',L/1000,filter_type,box_num,year),'-dpdf')
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
    
    ax = subplot(2,2,4);
    [~,h] = contourf(lon(opt_patch_lon), lat(opt_patch_lat),...
        nanmean(slhf_patch_bar+sshf_patch_bar,3)'-nanmean(A+B+C1+C2+C3+D+E1+E2,3)');
    title('$\overline{SSHF + SLHF}$ -$\Sigma A-E $ Wm$^{-2}$','interpreter','latex')
    format_fig(h,ax)
    
    figure(3)
    
    ax = subplot(3,4,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(dh_CTRL,3)');
    title('total:$\bar{\Delta h}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(dh_prime,3)');
    title('total:$\Delta h''$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(qo_patch(opt_patch_lon,opt_patch_lat,:),3)');
    title('total:$q_o$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,5);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(adh_CTRL,3)');
    title('total:$\bar{\alpha \Delta h}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,6);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(adh_prime,3)');
    title('total:$\alpha \Delta h''$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,7);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(qa_patch(opt_patch_lon,opt_patch_lat,:),3)');
    title('total:$q_a$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,8);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(qo_patch(opt_patch_lon,opt_patch_lat,:)-qa_patch(opt_patch_lon,opt_patch_lat,:),3)');
    title('total:$\Delta q$','interpreter','latex')
    format_fig(h,ax)
        
    ax = subplot(3,4,9);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(U_CTRL,3)');
    title('total:$\bar{U}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,10);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(U_prime,3)');
    title('total:$U''$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,11);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(To_prime,3)');
    title('total:$T_o''$','interpreter','latex')
    format_fig(h,ax)
    
     ax = subplot(3,4,12);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(DT_patch(opt_patch_lon,opt_patch_lat,:),3)');
    title('total:$\Delta T$','interpreter','latex')
    format_fig(h,ax)
    
    set(gcf,'position',[ 85          74        1231         731])
    update_figure_paper_size()
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/components_L_%d_%s_box%d_%d_synthetic',L/1000,filter_type,box_num,year),'-dpdf')
    
    figure(4)
    ax = subplot(2,1,1);
    [~,h] = contourf(lon_plot,lat_plot,(nanmean(To_prime.^2,3))','linestyle','none');
    c = get(gca,'clim');
hold on
[C,h] = contour(lon_plot,lat_plot,nanmean(SST,3)','w','linewidth',2);
    clabel(C,h,'color','w')
    set(gca,'clim',c)
    title('$(T_o'')^2$','interpreter','latex')
    set(gca,'ydir','normal','fontsize',15)
    colorbar
    xlabel('deg')
    ylabel('deg')
    
    ax = subplot(2,1,2);
    [~,h] = contourf(lon_plot,lat_plot,(nanmean(qs_prime.^2,3))','linestyle','none');
    c = get(gca,'clim');
    hold on
    [C,h] = contour(lon_plot,lat_plot,nanmean(SST,3)','w','linewidth',2);
    clabel(C,h,'color','w')
    set(gca,'clim',c)
    title('$(q_o'')^2$','interpreter','latex')
    set(gca,'ydir','normal','fontsize',15)
    colorbar
    xlabel('deg')
    ylabel('deg')
    set(gcf,'color','w','position',[61 221 1359 581])
    
    pause(1)
    drawnow
    pause(1)
    update_figure_paper_size()
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/To_L_%d_%s_box%d_%d_synthetic',L/1000,filter_type,box_num,year),'-dpdf')
    
    
    figure(5)
    ax = subplot(2,1,1);
    To2 = nanmean(To_prime.^2,3);
    [To2_CTRL,To2_prime] = FFT2D_filter(To2,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,To2_prime','linestyle','none');
    c = get(gca,'clim');
    hold on
    [C,h] = contour(lon_plot,lat_plot,nanmean(SST,3)','w','linewidth',2);
    clabel(C,h,'color','w')
    set(gca,'clim',c)
    title('$\overline{(T_o'')^2}$','interpreter','latex')
    set(gca,'ydir','normal','fontsize',15)
    colorbar
    xlabel('deg')
    ylabel('deg')
    
    ax = subplot(2,1,2);
    qo2 = nanmean(qs_prime.^2,3);
    [qo2_CTRL,qo2_prime] = FFT2D_filter(qo2,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,qo2_prime','linestyle','none');
    c = get(gca,'clim');
    hold on
    [C,h] = contour(lon_plot,lat_plot,nanmean(SST,3)','w','linewidth',2);
    clabel(C,h,'color','w')
    set(gca,'clim',c)
    title('$\overline{(q_o'')^2}$','interpreter','latex')
    set(gca,'ydir','normal','fontsize',15)
    colorbar
    xlabel('deg')
    ylabel('deg')
    set(gcf,'color','w','position',[61 221 1359 581])
    update_figure_paper_size()
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/To2_qo2_L_%d_%s_box%d_%d_synthetic',L/1000,filter_type,box_num,year),'-dpdf')
    
  figure(6)
     ax = subplot(3,2,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(A,3)');
    title('A: time mean of product $\rho_a C_D\overline{\overline{U}\,\overline{\Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,2,3);
    approx_A = rho_a*CD*nanmean(U_CTRL,3).*nanmean(dh_CTRL,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_A');
    title('approx A: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,2,5);
    [approx_A_LP,~] = FFT2D_filter(approx_A,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_A_LP');
    title('approx A: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,2,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(B,3)');
    title('B:$\rho_a C_D\overline{U''\Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,2,4);
    approx_B = rho_a*CD*nanmean(U_prime,3).*nanmean(dh_prime,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_B');
    title('approx B: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,2,6);
    [approx_B_LP,~] = FFT2D_filter(approx_B,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_B_LP');
    title('approx B: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    figure(7)
    ax = subplot(3,3,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C1,3)');
    title('C1: time mean of product $\rho_a C_D\overline{\overline{U} T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,4);
    approx_C1 = rho_a*CD*nanmean(U_CTRL,3).*nanmean(To_prime,3).*nanmean(adh_prime,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(approx_C1,3)');
    title('C1: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,7);
    [approx_C1_LP,~] = FFT2D_filter(approx_C1,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C1,3)');
    title('C1: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    
    ax = subplot(3,3,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C2,3)');
    title('C2: time mean of product $\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,5);
    approx_C2 = rho_a*CD*nanmean(U_prime,3).*nanmean(To_prime,3).*nanmean(adh_CTRL,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_C2');
    title('C2: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,8);
    [approx_C2_LP,~] = FFT2D_filter(approx_C2,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_C2_LP');
    title('C2: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(C3,3)');
    title('C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,6);
    approx_C3 = rho_a*CD*nanmean(U_prime,3).*nanmean(To_prime,3).*nanmean(adh_prime,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_C3');
    title('C3: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,9);
    [approx_C3_LP,~] = FFT2D_filter(approx_C3,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_C3_LP');
    title('C3: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    

    figure(8)
    ax = subplot(3,3,1);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(D,3)');
    title('D:$\rho_a C_D\overline{U''\overline{\Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,4);
    approx_D = rho_a*CD*nanmean(U_prime,3).*nanmean(dh_CTRL,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_D');
    title('D: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,7);
    [approx_D_LP,~] = FFT2D_filter(approx_D,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_D_LP');
    title('D: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,2);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(E1,3)');
    title('E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,5);
    approx_E1 = rho_a*CD*nanmean(U_CTRL,3).*nanmean(To_prime,3).*nanmean(adh_CTRL,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_E1');
    title('E1: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,8);
    [approx_E1_LP,~] = FFT2D_filter(approx_E1,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_E1_LP');
    title('E1: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,3);
    [~,h] = contourf(lon_plot,lat_plot,nanmean(E2,3)');
    title('E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,6);
    approx_E2 = rho_a*CD*nanmean(U_CTRL,3).*nanmean(dh_prime,3);
    [~,h] = contourf(lon_plot,lat_plot,approx_E2');
    title('E2: product of time mean','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,3,9);
    [approx_E2_LP,~] = FFT2D_filter(approx_E2,dx,cf,debug_flag,lat_plot,lon_plot);
    [~,h] = contourf(lon_plot,lat_plot,approx_E2_LP');
    title('E2: low pass of product of time mean','interpreter','latex')
    format_fig(h,ax)
    

    pause(1)
    drawnow
    pause(1)

    
    
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
