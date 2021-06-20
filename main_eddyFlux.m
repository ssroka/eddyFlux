clear;close all;clc

% move to this file's location
fileloc = mfilename('fullpath');
filedir = fileloc(1:max(strfind(fileloc,'/')));
cd(filedir)

% add paths to various filter locations
addpath('~/Documents/MATLAB/util/')
addpath('~/MIT/Research/eddyFlux/filter/')
addpath('~/MIT/Research/eddyFlux/filter/fft_filter/')
addpath('~/MIT/Research/eddyFlux/')
addpath('~/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

% add path to ERA5 data
addpath('~/MIT/Research/eddyFlux/ERA5_data/')
data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

%% USER input
year_vec = [2003 2007];
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
box_num = 5;
intvl = 1; % look at every intvl'th ti mpepoint
model_str = 'alpha'; % 'alpha' 'beta' 'alphabeta'
load('env_const.mat'); % load rho_a and c_p_air

debug_flag = false;
% plot_flag = false;
calc_alpha_beta_CD_flag = true;
plot_model_ERA_err_flag = true;
abCD_factor = 1;
alpha_pos_flag = false;
beta_pos_flag = false;
fft_first_flag = true;
plot_all_alpha_beta_CD = false;
plot_all_ABC_vs_year = false;
plot_all_QandC = false;
plot_ABC_comp = false;
%% SETUP
if plot_all_alpha_beta_CD || plot_all_ABC_vs_year
    model_str_cell = {'alpha','beta','alphabeta'};
else
    model_str_cell = {model_str};
end

for ii_model_str = 1:length(model_str_cell)
    model_str = model_str_cell{ii_model_str};
    switch model_str
        case 'alpha'
            abCD0 =   [  0.01    0.0014]';
        case 'beta'
            abCD0 =   [  0.3     0.0014]';
        case 'alphabeta'
            abCD0 =   [  0.01  0.3 0.0014]';
    end
    
    switch box_num
        case 1
            box_opt = [36 41.5; 143 152];
        case 2
            box_opt = [30 44.5; 148 169];
        case 3
            box_opt = [30 41.5; 142.5 169];
        case 4
            box_opt = [30 41.5; 142.5 153];
        case 5
            box_opt = [30 41.5; 153 169];
    end
    
    if alpha_pos_flag
        LB = [0;0;];
        UB = [Inf; Inf;];
        con_str = 'cons_';
    else
        LB = [-Inf; 0];
        UB = [Inf; Inf];
        con_str = '';
    end
    
    if beta_pos_flag
        LB = [0;0];
        UB = [Inf; Inf];
        con_str = 'cons_';
    else
        LB = [-Inf; 0];
        UB = [Inf; Inf];
        con_str = '';
    end
    
    if strcmp(filter_type,'fft')
        cf = (1/(2*L));
    end
    
    if fft_first_flag
        fft_str = '';
    else
        % fft_str = 'noFFT_'; USE calc_ABC_time_mean_then_fft if you want this
        % option
    end
    
    %% RUN
    setup_lat_lon_vec;
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
        'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');
    % only used to get the fft filter sampling rate;
    d_lat = abs(lat(2)-lat(1));
    m_per_deg = 111320;
    dx = abs(d_lat*m_per_deg);
    
    if ~(plot_all_alpha_beta_CD || plot_all_ABC_vs_year || plot_all_QandC || plot_ABC_comp)
        make_multipliers_eddyFlux(year_vec,L,data_src,filter_type,box_num,box_opt,cf,dx,debug_flag,model_str);
        get_CD_a_b_eddyFlux;
        calc_ABC_eddyFlux;
        plot_ABC_eddyFlux;
    else
        if plot_all_alpha_beta_CD
            plot_alpha_beta_CD_per_param
        end
        if plot_all_ABC_vs_year
            plot_ABC_terms_all_years
        end
        if plot_all_QandC
            plot_Qiu_Chen_03_19;
        end
        if plot_ABC_comp
            plot_ABC_comp_eddyFlux;
        end
    end
    
    
end

