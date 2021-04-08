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

fft_first_flag = true;

debug_flag = false;
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
    'lat','lon','patch_lat','patch_lon');

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

for i = 1:length(year_vec)
    year = year_vec(i);
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile)
    
    SST_prime = zeros(sum(opt_patch_lon),sum(opt_patch_lat),length(time));
    
    % coefficients
    load(sprintf('get_CD_alpha/opt_aCD_%sfilt_%s_L_%d_%d_to_%d',con_str,filter_type,L/1000,year_vec(1),year_vec(i)),'X');
    
    alpha_CD = X{i};
    
    as(i) = alpha_CD(1);
    aL(i) = alpha_CD(2);
    a(i) = mean([as aL]);
    
    CD_s(i) = alpha_CD(3);
    CD_L(i) = alpha_CD(4);
    CD(i) = mean([CD_s CD_L]);
    
end
%%
mean_CD_all_years = mean(CD);

figure(1);clf
yyaxis left
plot(year_vec,CD_s,'k','displayname','$C_D^s$','linewidth',2)
hold on
plot(year_vec,CD_L,'k--','displayname','$C_D^L$','linewidth',2)
set(gca,'ycolor','k')
ylabel('C_D')
xlabel('year')
yyaxis right
plot(year_vec,abs(CD_s-CD_L)./CD,'r','displayname','$\frac{|C_D^s-C_D^L|}{0.5(CD_s*CD_L)}$ each year','linewidth',2)
hold on
plot(year_vec,abs(mean_CD_all_years-CD_s)./mean_CD_all_years,'r--','displayname','$\frac{|CD_{avg}-C_D^s|}{CD_{avg}}$ ','linewidth',2)
plot(year_vec,abs(mean_CD_all_years-CD_L)./mean_CD_all_years,'r:','displayname','$\frac{|CD_{avg}-C_D^L|}{CD_{avg}}$ ','linewidth',2)


set(gca,'ycolor','r','fontsize',25)
legend('location','northoutside','interpreter','latex')
ylabel('relative error')
xlabel('year')
set(gcf,'color','w','position',[61 24 1359 778],'NumberTitle','off','Name',num2str(year))
update_figure_paper_size()
print(sprintf('imgs/CD_variability',L/1000,filter_type,year),'-dpdf')

%%
figure(2)
yyaxis left
plot(year_vec,as,'k','displayname','$\alpha_s$','linewidth',2)
hold on
plot(year_vec,aL,'k--','displayname','$\alpha_L$','linewidth',2)
legend('location','best','interpreter','latex')
yyaxis right
plot(year_vec,(as-aL)./a,'r','displayname','$\frac{CD_s-CD_L}{0.5(CD_s*CD_L)}$','linewidth',2)






