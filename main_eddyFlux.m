clear;close all;clc

% move to this file's location
fileloc = mfilename('fullpath');
filedir = fileloc(1:max(strfind(fileloc,'/')));
cd(filedir)

% data_base = '/Users/ssroka/MIT/Research/eddyFlux/';
data_base = '/Volumes/SydneySroka_Remy/eddyFlux/';

% add paths to various filter locations
addpath('~/Documents/MATLAB/util/')
addpath([data_base 'filter/'])
addpath([data_base 'filter/fft_filter/'])
addpath(data_base)
addpath([data_base 'get_CD_alpha'])
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')



%% USER input
year_vec = [2003];
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
box_num = 3;
er_box_num = 3;

reanalysis_src = 'NCEP';% 'ERA5' 'NCEP'

if strcmp(reanalysis_src,'ERA5')
    % add path to ERA5 data
    addpath([data_base 'ERA5_data/'])
    data_src =[data_base  'ERA5_data/'];
else
    % add path to NCEP data
    addpath([data_base 'NCEP_data/'])
    data_src =[data_base  'NCEP_data/'];
end

for abCD_fac_vec = [-10 -8 -5 -4 -2 -0.5 1 0.5 2 4 5 8 10]
intvl = 1; % look at every intvl'th ti mpepoint
model_str = 'beta'; % 'alpha' 'beta' 'alphabeta'
load('env_const.mat'); % load rho_a and c_p_air

debug_flag = false;
% plot_flag = false;
calc_alpha_beta_CD_flag = false;
plot_model_ERA_err_flag = false;
calc_model_ERA_rms_err_flag = true;
alpha_pos_flag = false;
beta_pos_flag = false;
fft_first_flag = true; % don't touch!!
plot_all_alpha_beta_CD = false;
plot_all_ABC_vs_year = false;
plot_all_QandC = false;
plot_ABC_comp = false;
%% SETUP
if strcmp(model_str,'alphabeta')
    abCD_factor = [abCD_fac_vec abCD_fac_vec 1]';
else
    abCD_factor = [abCD_fac_vec 1]';
end


if plot_all_alpha_beta_CD || plot_all_ABC_vs_year
%     model_str_cell = {'alpha','beta','alphabeta'};
    model_str_cell = {'beta','alphabeta'};
else
    model_str_cell = {model_str};
end

for ii_model_str = 1:length(model_str_cell)
    model_str = model_str_cell{ii_model_str};
    switch model_str
        case 'alpha'
            abCD0 =   [  0.01    0.0014]';
            model_str_title = '\alpha';
        case 'beta'
            abCD0 =   [  0.3     0.0014]';
            model_str_title = '\beta';
        case 'alphabeta'
            abCD0 =   [  0.01  0.3 0.0014]';
            model_str_title = '\alpha\beta';
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
    
    switch er_box_num
        case 3
            er_box = [30 41.5; 142.5 169];% [35.75 41.5; 142.5 155.7500];
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
    if all(abs(abCD_factor-1)<1e-10)
        abCD_fac_str = '';
    else
        abCD_fac_str = ['_'];
        for abCD_ii = 1:length(abCD_factor)
            abCD_fac_str = [abCD_fac_str strrep(num2str(abCD_factor(abCD_ii)),'.','_') '_'];
        end
        abCD_fac_str = abCD_fac_str(1:end-1);
    end
    %% RUN
    setup_lat_lon_vec;
    if strcmp(reanalysis_src,'ERA5')
        load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
            'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');
    else
        load(sprintf('%sNCEP_patch_data_%d.mat',data_src,year_vec(1)),...
            'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');
    end
    % only used to get the fft filter sampling rate;
    d_lat = abs(lat(2)-lat(1));
    m_per_deg = 111320;
    dx = abs(d_lat*m_per_deg);
    
    if ~(plot_all_alpha_beta_CD || plot_all_ABC_vs_year || plot_all_QandC || plot_ABC_comp || calc_model_ERA_rms_err_flag)
        make_multipliers_eddyFlux(year_vec,L,data_src,filter_type,box_num,box_opt,cf,dx,debug_flag,model_str,er_box,reanalysis_src);
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
        if calc_model_ERA_rms_err_flag
            calc_rms_error_model_ERA5;
        end
    end
    
    
end
end

