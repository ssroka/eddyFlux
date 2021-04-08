clear;clc;close all
year1 = 2003;
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

fft_first_flag = true;

param_3_str = true;

debug_flag = false;

box_num = 1;

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
    %     box_opt = [30 42; 144 168];
%     box_opt = [30 44.5; 148 169];
%     box_opt = [36 41.5; 143 152];
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
    'lat','lon','patch_lat','patch_lon');

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

if param_3_str
    param_num_str = '_3param';
else
    param_num_str = '';

end

for i = 1:length(year_vec)
    year = year_vec(i);
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile)
    
    SST_prime = zeros(sum(opt_patch_lon),sum(opt_patch_lat),length(time));
    
    % coefficients
    if param_3_str
        load(sprintf('get_CD_alpha/opt_aCD_%sfilt_%s_L_%d%s_box%d_%d',con_str,filter_type,L/1000,param_num_str,box_num,year_vec(i)),'X');

        alpha_CD = X{i};
        
        as = alpha_CD(1);
        aL = alpha_CD(2);
        
        CD = alpha_CD(3);
    else
        load(sprintf('get_CD_alpha/opt_aCD_%sfilt_%s_L_%d%s_box%d_%d',con_str,filter_type,L/1000,param_num_str,box_num,year_vec(i)),'X');

        alpha_CD = X{i};
        
        as = alpha_CD(1);
        aL = alpha_CD(2);
        
        CD_s = alpha_CD(3);
        CD_L = alpha_CD(4);
        CD = mean([CD_s CD_L]);
        
    end

    
    p = size(SST_patch,3); % re calculate for leap year
    
    ZEROS = zeros(m,n,p);
    
    
    A = ZEROS;
    B = ZEROS;
    C1 = ZEROS;
    C2 = ZEROS;
    C3 = ZEROS;
    D = ZEROS;
    E1 = ZEROS;
    E2 = ZEROS;
    sshf_patch_bar = ZEROS;
    slhf_patch_bar = ZEROS;
    adh_CTRL = ZEROS;
    adh_prime = ZEROS;
    U_CTRL = ZEROS;
    U_prime = ZEROS;
    To_prime = ZEROS;
    dh_CTRL =  ZEROS;
    dh_prime =  ZEROS;
    qs_prime =  ZEROS;
    SST =  ZEROS;
    
    count = 1;
    fprintf('\n')
    for tt = 1:intvl:p % time points
        fprintf(' processing snapshot %d of %d\n',tt,p)
        
        if strcmp(filter_type,'boxcar')
            [SST_patch_CTRL,SST_prime] = boxcar_filter(SST_patch(:,:,tt),M);
            %             [P0_patch_CTRL,~] = boxcar_filter(P0_patch(:,:,tt),M);
            [q_diff_CTRL,~] = boxcar_filter(qo_patch(:,:,tt)-qa_patch(:,:,tt),M);
            [DT_patch_CTRL,~] = boxcar_filter(DT_patch(:,:,tt),M);
            [U_mag_CTRL,~] = boxcar_filter(U_mag(:,:,tt),M);
            U_mag_prime = U_mag(:,:,tt) - U_mag_CTRL;
            
        elseif strcmp(filter_type,'fft')
%             debug_flag = false;

            % SST
            [SST_patch_CTRL,SST_prime] = FFT2D_filter(SST_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            
            % qo
            [qo_CTRL,qo_prime] = FFT2D_filter(qo_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            %             SST_prime = SST_patch(opt_patch_lon,opt_patch_lat,tt) - SST_patch_CTRL;
            
            Lv = SW_LatentHeat(SST_patch(opt_patch_lon,opt_patch_lat,tt),'K',salinity,'ppt');
            
            % h
            h_diff = Lv.*(qo_patch(opt_patch_lon,opt_patch_lat,tt)-qa_patch(opt_patch_lon,opt_patch_lat,tt)) +...
                c_p_air.*(DT_patch(opt_patch_lon,opt_patch_lat,tt));
            [h_diff_CTRL,h_diff_prime] = FFT2D_filter(h_diff,dx,cf,debug_flag,lat_plot,lon_plot);
            %             h_diff_prime = h_diff - h_diff_CTRL;
            
            % alpha h
            ah_diff = aL.*Lv.*(qo_patch(opt_patch_lon,opt_patch_lat,tt)-qa_patch(opt_patch_lon,opt_patch_lat,tt)) +...
                as.*c_p_air.*(DT_patch(opt_patch_lon,opt_patch_lat,tt));
            [ah_diff_CTRL,ah_diff_prime] = FFT2D_filter(ah_diff,dx,cf,debug_flag,lat_plot,lon_plot);
            %             ah_diff_prime = ah_diff - ah_diff_CTRL;
            
            % U
            [U_mag_CTRL,U_mag_prime] = FFT2D_filter(U_mag(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            %             U_mag_prime = U_mag(opt_patch_lon,opt_patch_lat,tt) - U_mag_CTRL;
            
            A_full = rho_a*CD*U_mag_CTRL.*h_diff_CTRL;
            B_full = rho_a*CD*U_mag_prime.*h_diff_prime;
            C1_full = rho_a*CD*U_mag_CTRL.*SST_prime.*ah_diff_prime;
            C2_full = rho_a*CD*U_mag_prime.*SST_prime.*ah_diff_CTRL;
            C3_full = rho_a*CD*U_mag_prime.*SST_prime.*ah_diff_prime;
            D_full = rho_a*CD*U_mag_prime.*h_diff_CTRL;
            E1_full = rho_a*CD*SST_prime.*U_mag_CTRL.*ah_diff_CTRL;
            E2_full = rho_a*CD*U_mag_CTRL.*h_diff_prime;
            
            
            if fft_first_flag
%                 debug_flag = 1;
                A(:,:,count) = FFT2D_filter(A_full,dx,cf,debug_flag,lat_plot,lon_plot);
                B(:,:,count) = FFT2D_filter(B_full,dx,cf,debug_flag,lat_plot,lon_plot);
                C1(:,:,count) = FFT2D_filter(C1_full,dx,cf,debug_flag,lat_plot,lon_plot);
                C2(:,:,count) = FFT2D_filter(C2_full,dx,cf,debug_flag,lat_plot,lon_plot);
                C3(:,:,count) = FFT2D_filter(C3_full,dx,cf,debug_flag,lat_plot,lon_plot);
                D(:,:,count) = FFT2D_filter(D_full,dx,cf,debug_flag,lat_plot,lon_plot);
                E1(:,:,count) = FFT2D_filter(E1_full,dx,cf,debug_flag,lat_plot,lon_plot);
                E2(:,:,count) = FFT2D_filter(E2_full,dx,cf,debug_flag,lat_plot,lon_plot);
                sshf_patch_bar(:,:,count) = FFT2D_filter(sshf_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
                slhf_patch_bar(:,:,count) = FFT2D_filter(slhf_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
                fft_str = '';
                
            else
                A(:,:,count) = A_full;
                B(:,:,count) = B_full;
                C1(:,:,count) = C1_full;
                C2(:,:,count) = C2_full;
                C3(:,:,count) = C3_full;
                D(:,:,count) = D_full;
                E1(:,:,count) = E1_full;
                E2(:,:,count) = E2_full;
                sshf_patch_bar(:,:,count) = sshf_patch(opt_patch_lon,opt_patch_lat,tt);
                slhf_patch_bar(:,:,count) = slhf_patch(opt_patch_lon,opt_patch_lat,tt);
                
                fft_str = 'noFFT_';
            end
            
            
            %             if tt == 1
            %                 study_ABC
            %             end
            
            dh_CTRL(:,:,count) = h_diff_CTRL;
            dh_prime(:,:,count) = h_diff_prime;
            adh_CTRL(:,:,count) = ah_diff_CTRL;
            adh_prime(:,:,count) = ah_diff_prime;
            U_CTRL(:,:,count) = U_mag_CTRL;
            U_prime(:,:,count) = U_mag_prime;
            To_prime(:,:,count) = SST_prime;
            qs_prime(:,:,count) = qo_prime;
            SST(:,:,count) = SST_patch(opt_patch_lon,opt_patch_lat,tt);
            
        elseif strcmp(filter_type(1:7),'lanczos')
            [SST_patch_CTRL,SST_prime] = lanczos_filter(SST_patch(:,:,tt),dx,cf);
            %             [P0_patch_CTRL,~] = lanczos_filter(P0_patch(:,:,tt),dx,cf);
            [q_diff_CTRL,~] = lanczos_filter(qo_patch(:,:,tt)-qa_patch(:,:,tt),dx,cf);
            [DT_patch_CTRL,~] = lanczos_filter(DT_patch(:,:,tt),dx,cf);
            [U_mag_CTRL,~] = lanczos_filter(U_mag(:,:,tt),dx,cf);
            U_mag_prime = U_mag(:,:,tt) - U_mag_CTRL;
            
        end
        
        
        
        count = count + 1;
    end
    save(sprintf('ABC_terms_%d_%sfilt_%s%s%s_box%d_%d',L/1000,con_str,fft_str,filter_type,param_num_str,box_num,year),...
        'A','B','C1','C2','C3','D','E1','E2','slhf_patch_bar','sshf_patch_bar','SST',...
        'dh_CTRL','dh_prime','adh_CTRL','adh_prime','U_CTRL','U_prime','To_prime','qs_prime',...
        'dx','cf','debug_flag','lat_plot','lon_plot','as','aL','CD')
    fprintf('\nsaving A,B and C for %d\n',year)
    
end


% get_eddy_contribution_plot;
