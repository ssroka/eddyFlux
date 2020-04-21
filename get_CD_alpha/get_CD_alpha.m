
clear;close all;clc
addpath('~/Documents/MATLAB/util/')
addpath('~/MIT/Research/eddyFlux/filter/')
addpath('~/MIT/Research/eddyFlux/ERA5_data/')

calc_CD_alpha_flag = false;
plot_flag = true;

%     % Kurishio
%     lat_bnds = [25 45];
%     lon_bnds = [130 170];

err_box_lat = [32 38];
err_box_lon = [140 160];

L = 500000; % m

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

year_vec = 2003:2007;

% 2003
aCd0(:,1) = [  0.0136; -0.0447; 0.0014; 0.0023];

% 2004
aCd0(:,2) = [ 0.0180; -0.0501; 0.0014; 0.0024];

% 2005
aCd0(:,3) = [ 0.0167; -0.0458; 0.0014; 0.0023];

% 2006
aCd0(:,4) = [ 0.0054;-0.0507; 0.0014; 0.0022];

% 2007
aCd0(:,5) = [ 0.0133; -0.0481; 0.0014; 0.0024];


%% parameters for optimization

LB = [-Inf; -Inf; 0; 0];
UB = [Inf; Inf; Inf; Inf];

options = optimoptions(@fmincon,'Display','iter','steptolerance',1e-5);

%% begin optimizing for CD and alpha

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'lat','lon','patch_lat','patch_lon');

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));


X = cell(length(year_vec),1);
FFINAL = cell(length(year_vec),1);

for i = 1:length(year_vec)
    
    if calc_CD_alpha_flag
        
        [x,ffinal] = fmincon(@(aCd) optimize_alpha_CD(aCd,year_vec(i),err_box_bnds_lat,err_box_bnds_lon),aCd0(:,year_vec(i)-2002),[],[],[],[],LB,UB,[],options);
        
        X{i} = x;
        FFINAL{i}  =ffinal;
        save(sprintf('opt_alpha_CD_%d_to_%d',year_vec(1),year_vec(i)),'X','FFINAL');
        fprintf('saved %d\n',year_vec(i))
        [aCd0(:,year_vec(i)-2002) X{i}]
        
    elseif plot_flag
        
        year = year_vec(i);
        
        filename = sprintf('Qs_QL_optimization_data_%d',year_vec(i));
        load(filename,'as_multiplier','aL_multiplier','SST_prime','sshf_patch','slhf_patch')

        
        lat_er = lat(patch_lat);
        lat_er = lat_er(err_box_bnds_lat);
        lon_er = lon(patch_lon);
        lon_er = lon_er(err_box_bnds_lon);
        
        t_range = 1:size(sshf_patch,3);
        
        load(sprintf('opt_alpha_CD_%d_to_%d',year_vec(1),year_vec(end)),'X');
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








































