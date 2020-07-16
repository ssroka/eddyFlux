
clear;close all;clc
addpath('~/Documents/MATLAB/util/')
addpath('~/MIT/Research/eddyFlux/filter/')
addpath('~/MIT/Research/eddyFlux/ERA5_data/')

calc_CD_alpha_flag = true;
plot_flag = false;

alpha_pos_flag = false;

%     % Kurishio
%     lat_bnds = [25 45];
%     lon_bnds = [130 170];


L = 500000; % m

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'

year_vec = [2003 2007];

% % 2003
% aCd0(:,1) = [  0.0105; 0.0165; 0.0014; 0.0015];
% 
% % 2004
% aCd0(:,2) = [ 0.0114; 0.0145; 0.0014; 0.0015];
% 
% % 2005
% aCd0(:,3) = [ 0.0106; 0.0126; 0.0014; 0.0015];
% 
% % 2006
% aCd0(:,4) = [ -0.0004; 0.0056; 0.0014; 0.0015];
% 
% % 2007
% aCd0(:,5) = [ 0.0125; 0.0115; 0.0014; 0.0015];

aCd0 = [0.0038 0.0045 0.0014 0.0015]'; % lanczos
% aCd0 = [  0.0105; 0.0165; 0.0014; 0.0015]; % boxcar
aCd0 = [0.0038 0.0045 0.0014 0.0015]'; % fft


%% parameters for optimization

if strcmp(filter_type,'boxcar')
    box_limits = [32 38; 140 160];
elseif strcmp(filter_type,'fft')
    cf = (1/(2*L));
    box_limits = [30 42; 148 168];
elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
    box_limits = [32 38; 140 160];
end

if alpha_pos_flag
    LB = [0;0; 0; 0];
    UB = [Inf; Inf; Inf; Inf];
    con_str = 'cons_';
else
    LB = [-Inf; -Inf; 0; 0];
    UB = [Inf; Inf; Inf; Inf];
    con_str = '';
end

options = optimoptions(@fmincon,'Display','iter','steptolerance',1e-5);

%% begin optimizing for CD and alpha

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'lat','lon','patch_lat','patch_lon');

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));

box_lat = lat>=box_limits(1,1) & lat<=box_limits(1,2);
box_lon = lon>=box_limits(2,1) & lon<=box_limits(2,2);

X = cell(length(year_vec),1);
FFINAL = cell(length(year_vec),1);

for i = 1:length(year_vec)
    
    if calc_CD_alpha_flag
        
        [x,ffinal] = fmincon(@(aCd) optimize_alpha_CD(aCd,year_vec(i),err_box_bnds_lat,err_box_bnds_lon,L,filter_type),aCd0,[],[],[],[],LB,UB,[],options);
        
        X{i} = x;
        FFINAL{i}  =ffinal;
        save(sprintf('opt_aCD_%sfilt_%s_L_%d_%d_to_%d',con_str,filter_type,L/1000,year_vec(1),year_vec(i)),'X','FFINAL');
        fprintf('saved %d\n',year_vec(i))
        [aCd0 X{i}]
        
    end
    if plot_flag
        
        year = year_vec(i);
        filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_%d',L/1000,filter_type,year_vec(i));
%         filename = sprintf('Qs_QL_optimization_data_%d',year_vec(i));
        load(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')

        
        lat_er = lat(patch_lat);
        lat_er = lat_er(err_box_bnds_lat);
        lon_er = lon(patch_lon);
        lon_er = lon_er(err_box_bnds_lon);
        
        t_range = 1:size(sshf_patch,3);
        
        load(sprintf('opt_aCD_%sfilt_%s_L_%d_%d_to_%d',con_str,filter_type,L/1000,year_vec(1),year_vec(i)),'X','FFINAL');

        alpha_CD = X{i};
        
        as = alpha_CD(1);
        aL = alpha_CD(2);
        
        CD_s = alpha_CD(3);
        CD_L = alpha_CD(4);
        
        sshf_model = CD_s.*(1+as.*SST_prime).*as_multiplier;
        slhf_model = CD_L.*(1+aL.*SST_prime).*aL_multiplier;
        
        figure(year)
        cmp_ERA5_model;
        
    end
    
    
end






