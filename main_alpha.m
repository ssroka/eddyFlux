
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

% USER input
year_vec = [2003];
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
box_num = 2;





debug_flag = false;
calc_CD_alpha_flag = false;
plot_flag = true;
alpha_pos_flag = false;

% initial guess for alpha and CD
aCd0 =   [  0.0104  0.0135    0.0014]'; % fft L = 250 km 

param_num_str = '_3param';

%%

% which box to use
switch box_num
    case 1
        box_opt = [36 41.5; 143 152];
    case 2
        box_opt = [30 44.5; 148 169];
end

if strcmp(filter_type,'boxcar')
    box_opt = [25 45; 130 170];
elseif strcmp(filter_type,'fft')
    cf = (1/(2*L));
elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
    box_opt = [25 45; 130 170];
end

if alpha_pos_flag
    LB = [0;0; 0];
    UB = [Inf; Inf; Inf];
    con_str = 'cons_';
else
    LB = [-Inf; -Inf; 0];
    UB = [Inf; Inf; Inf];
    con_str = '';
end


%%
cd get_CD_alpha
% make_multipliers(year_vec,L,data_src,filter_type,box_num,box_opt,cf,debug_flag);
% % make_multipliers
get_CD_alpha_3param;

cd ..
calc_ABC
plot_ABC













