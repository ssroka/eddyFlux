clear;clc;close all

year_vec = [2003:2019];

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


param_3_str = true;

fft_first_flag = true; % Do nuot change, USE calc_ABC_time_mean_then_fft if you want this
% option to be false

add_SSH_var = true;
%% filter set up

load('env_const.mat')

files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),'SST_patch','lat','lon','patch_lat','patch_lon');


files_for_size.lat = double(files_for_size.lat);
files_for_size.lon = double(files_for_size.lon);

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
    %     box_opt = [30 42; 144 168];
    box_opt = [30 44.5; 148 169];
elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
end

% this is the same for every year
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');


lat = double(lat);
lon = double(lon);

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

totA = (abs(lat_plot(end)-lat_plot(1))*m_per_deg)*(abs(lon_plot(end)-lon_plot(1))*m_per_deg);

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

ZERO = zeros(length(year_vec),1);
avg_A  = ZERO;
avg_B  = ZERO;
avg_C1 = ZERO;
avg_C2 = ZERO;
avg_C3 = ZERO;
avg_D  = ZERO;
avg_E1 = ZERO;
avg_E2 = ZERO;

for i = 1:length(year_vec)
    year = year_vec(i);
    
    load(sprintf('ABC_terms_%d_%sfilt_%s%s%s_%d',L/1000,con_str,fft_str,filter_type,param_num_str,year))
    
    avg_A(i)  = sum(sum(nanmean(A,3)'*dx*dy))/totA;
    avg_B(i)  = sum(sum(nanmean(B,3)'*dx*dy))/totA;
    avg_C1(i) = sum(sum(nanmean(C1,3)'*dx*dy))/totA;
    avg_C2(i) = sum(sum(nanmean(C2,3)'*dx*dy))/totA;
    avg_C3(i) = sum(sum(nanmean(C3,3)'*dx*dy))/totA;
    avg_D(i)  = sum(sum(nanmean(D,3)'*dx*dy))/totA;
    avg_E1(i) = sum(sum(nanmean(E1,3)'*dx*dy))/totA;
    avg_E2(i) = sum(sum(nanmean(E2,3)'*dx*dy))/totA;
    

    
end


    figure(1)
    ax = subplot(2,4,1);
    h = plot(year_vec,avg_A,'linewidth',2);
    title('A:$\rho_a C_D\overline{U}\overline{\Delta h}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,2);
    h = plot(year_vec,avg_B,'linewidth',2);
    title('B:$\rho_a C_D\overline{U''\Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,3);
    h = plot(year_vec,avg_C1,'linewidth',2);
    title('C1:$\rho_a C_D\overline{T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,4);
    h = plot(year_vec,avg_C2,'linewidth',2);
    title('C2:$\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,5);
    h = plot(year_vec,avg_C3,'linewidth',2);
    title('C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
     ax = subplot(2,4,6);
    h = plot(year_vec,avg_D,'linewidth',2);
    title('D:$\rho_a C_D\overline{U''\overline{\Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,7);
    h = plot(year_vec,avg_E1,'linewidth',2);
    title('E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,8);
    h = plot(year_vec,avg_E2,'linewidth',2);
    title('E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
        
    if add_SSH_var
        SSH_data = load(sprintf('var_SSH_%d_%d',year_vec(1),min(year_vec(end),2018)),'var_SSH','year_vec');
    for i = 1:8
        subplot(2,4,i)
        yyaxis right
        plot(SSH_data.year_vec,SSH_data.var_SSH,'r--','linewidth',2,'displayname','var(SSH)')
        set(gca,'ycolor','r')
    end
    
    
    cc(1) = corr(SSH_data.var_SSH,avg_A(1:end-1));
    cc(2) = corr(SSH_data.var_SSH,avg_B(1:end-1));
    cc(3) = corr(SSH_data.var_SSH,avg_C1(1:end-1));
    cc(4) = corr(SSH_data.var_SSH,avg_C2(1:end-1));
    cc(5) = corr(SSH_data.var_SSH,avg_C3(1:end-1));
    cc(6) = corr(SSH_data.var_SSH,avg_D(1:end-1));
    cc(7) = corr(SSH_data.var_SSH,avg_E1(1:end-1));
    cc(8) = corr(SSH_data.var_SSH,avg_E2(1:end-1));
    end
    
    set(gcf,'color','w','position',[1 215 1439 587],'NumberTitle','off','Name',num2str(year))
    update_figure_paper_size()
    print(sprintf('imgs/ABC_L_%d%s_%s_%d_%d',L/1000,param_num_str,filter_type,year_vec(1),year_vec(end)),'-dpdf')
 


function [] = format_fig(h,plt_num,max_val,min_val)
set(gca,'fontsize',15)
xlabel('year')
ylabel('$W/m^2$','interpreter','latex')
set(h,'linewidth',2,'color','k')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end


% get_eddy_contribution_plot;
