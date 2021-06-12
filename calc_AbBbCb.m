% clear;clc;close all
% year1 = 2003;
% year_vec = [2003];
% 
% data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';
% 
% addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
% addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/lanczosfilter/')
% addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/fft_filter/')
% addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
% addpath('/Users/ssroka/Documents/MATLAB/util/')
% addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')
% 
% L = 250000; % m
% 
% intvl = 1; % look at every intvl'th timpepoint
% 
% beta_pos_flag = false;
% 
% filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
% 
% fft_first_flag = true;
% 
% debug_flag = false;
% 
% box_num = 1;

%% filter set up

load('env_const.mat')

files_for_size =  load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),'SST_patch','lat','lon','patch_lat','patch_lon');

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

% switch box_num
%     case 1
%         box_opt = [36 41.5; 143 152];
%     case 2
%         box_opt = [30 44.5; 148 169];
% end

% this is the same for every year
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon');

file_for_box = load(sprintf('beta_Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year_vec(1)),'box_opt');

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

if beta_pos_flag
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
        load(sprintf('opt_bCD_%sfilt_%s_L_%d_box%d_%d',con_str,filter_type,L/1000,box_num,year_vec(i)),'X');

        beta_CD = X{i};
        

        b = beta_CD(2);
        CD = beta_CD(1);


    p = size(SST_patch,3); % re calculate for leap year
    
    ZEROS = zeros(m,n,p);
    
    
    Ab = ZEROS;
    Bb = ZEROS;
    Cb = ZEROS;
    Db = ZEROS;
    sshf_patch_bar = ZEROS;
    slhf_patch_bar = ZEROS;
    U_bar = ZEROS;
    U_prime = ZEROS;
    To_prime = ZEROS;
    dh_bar =  ZEROS;
    dh_prime =  ZEROS;
    Dq_prime =  ZEROS;
    SST =  ZEROS;
    
    
    filename = sprintf('beta_Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year);
    load(filename,'bs_multiplier','bL_multiplier','SST_prime','U_bar','sshf_patch','slhf_patch')


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
            [SST_bar_tt,SST_prime_tt] = FFT2D_filter(SST_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            
            % DT
            [DT_bar_tt,DT_prime_tt] = FFT2D_filter(DT_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            
            % qo
            [qo_CTRL_tt,qo_prime_tt] = FFT2D_filter(qo_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            
            % Dq
            Dq = (qo_patch(opt_patch_lon,opt_patch_lat,tt)-qa_patch(opt_patch_lon,opt_patch_lat,tt));
            [Dq_bar_tt,Dq_prime_tt] = FFT2D_filter(Dq,dx,cf,debug_flag,lat_plot,lon_plot);
            
            Lv = SW_LatentHeat(SST_patch(opt_patch_lon,opt_patch_lat,tt),'K',salinity,'ppt');
            
            % h
            CD_h_diff = CD*(bs_multiplier(:,:,tt) + bL_multiplier(:,:,tt));
%             CD_h_diff = rho_a*CD.*Lv.*Dq +...
%                 rho_a*CD.*c_p_air.*(DT_patch(opt_patch_lon,opt_patch_lat,tt));
            [CD_h_diff_bar_tt,CD_h_diff_prime_tt] = FFT2D_filter(CD_h_diff,dx,cf,debug_flag,lat_plot,lon_plot);
            %             h_diff_prime = h_diff - h_diff_CTRL;
            
          
            % U
            [U_bar_tt,U_prime_tt] = FFT2D_filter(U_mag(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
            %             U_mag_prime = U_mag(opt_patch_lon,opt_patch_lat,tt) - U_mag_CTRL;
            
            Ab_full = U_bar_tt.*CD_h_diff;
            Bb_full = U_bar_tt.*CD_h_diff_prime_tt;
            Cb_full = b.*SST_prime_tt.*CD_h_diff_prime_tt;
            Db_full = b.*SST_prime_tt.*CD_h_diff_bar_tt;

            
            if fft_first_flag
%                 debug_flag = 1;
                Ab(:,:,count) = FFT2D_filter(Ab_full,dx,cf,debug_flag,lat_plot,lon_plot);
                Bb(:,:,count) = FFT2D_filter(Bb_full,dx,cf,debug_flag,lat_plot,lon_plot);
                Cb(:,:,count) = FFT2D_filter(Cb_full,dx,cf,debug_flag,lat_plot,lon_plot);
                Db(:,:,count) = FFT2D_filter(Db_full,dx,cf,debug_flag,lat_plot,lon_plot);
                sshf_patch_bar(:,:,count) = FFT2D_filter(sshf_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
                slhf_patch_bar(:,:,count) = FFT2D_filter(slhf_patch(opt_patch_lon,opt_patch_lat,tt),dx,cf,debug_flag,lat_plot,lon_plot);
                fft_str = '';
                
            else
                A(:,:,count) = Ab_full;
                Bb(:,:,count) = Bb_full;
                Cb(:,:,count) = Cb_full;
                Db(:,:,count) = Db_full;

                sshf_patch_bar(:,:,count) = sshf_patch(opt_patch_lon,opt_patch_lat,tt);
                slhf_patch_bar(:,:,count) = slhf_patch(opt_patch_lon,opt_patch_lat,tt);
                
                fft_str = 'noFFT_';
            end
            
            
            %             if tt == 1
            %                 study_ABC
            %             end
            
            dh_bar(:,:,count) = CD_h_diff_bar_tt;
            dh_prime(:,:,count) = CD_h_diff_prime_tt;
            U_bar(:,:,count) = U_bar_tt;
            U_prime(:,:,count) = U_prime_tt;
            To_prime(:,:,count) = SST_prime_tt;
            Dq_prime(:,:,count) = Dq_prime_tt;
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
    save(sprintf('AbBbCb_terms_%d_%sfilt_%s%s_box%d_%d',L/1000,con_str,fft_str,filter_type,box_num,year),...
        'Ab','Bb','Cb','Db','slhf_patch_bar','sshf_patch_bar','SST',...
        'dh_bar','dh_prime','U_bar','U_prime','To_prime','Dq_prime',...
        'dx','cf','debug_flag','lat_plot','lon_plot','b','CD')
    fprintf('\nsaving Ab,Bb and Cb for %d\n',year)
     
end


% get_eddy_contribution_plot;
