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

% add path to ERA5 data
addpath([data_base 'ERA5_data/'])
data_src =[data_base  'ERA5_data/'];

%% USER input
year_vec = [2003:2018];
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
box_num = 3;
abCD_fac_vec = [1];
intvl = 1; % look at every intvl'th ti mpepoint
model_str = 'beta'; % 'alpha' 'beta' 'alphabeta'
load('env_const.mat'); % load rho_a and c_p_air
debug_flag = false;


%% SETUP
if strcmp(model_str,'alphabeta')
    abCD_factor = [abCD_fac_vec abCD_fac_vec 1]';
else
    abCD_factor = [abCD_fac_vec 1]';
end
if strcmp(filter_type,'fft')
    cf = (1/(2*L));
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
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');
% only used to get the fft filter sampling rate;
d_lat = abs(lat(2)-lat(1));
m_per_deg = 111320;
dx = abs(d_lat*m_per_deg);

SST_prime_total = zeros(488,length(year_vec));


for i = 1:length(year_vec)
    year = year_vec(i);
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile,'SST_patch','time')
    SST_box = SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:);
    SST_timemean = nanmean(SST_box,3)';
    
    [SST_time_mean_bar] = FFT2D_filter(SST_timemean,dx,cf,debug_flag,lon_box,lat_box);
    SST_time_mean_prime = SST_timemean - SST_time_mean_bar;
    
    filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year);
    load(filename,'SST_prime')
    

    
    if i == 1
        Toprime_stdev = zeros(length(time),length(year_vec));
        Toprime_mean = zeros(length(time),length(year_vec));
        SST_prime_total = zeros(size(SST_prime,1)*size(SST_prime,2)*488,length(year_vec));
        SST_prime_mean = zeros(size(SST_prime,1)*size(SST_prime,2),length(year_vec));
    end
    SST_prime_total(1:numel(SST_prime),i) = SST_prime(:);
    SST_prime_mean(1:size(SST_prime,1)*size(SST_prime,2),i) = reshape(nanmean(SST_prime,3),size(SST_prime,1)*size(SST_prime,2),1);

end
SST_prime_total_1vec = SST_prime_total(:);
SST_prime_total_1vec(SST_prime_total_1vec==0) = [];
%%
figure(1)
histogram(SST_prime_mean(:),50,'normalization','pdf','displayname','$\langle T_o'' \rangle $')
hold on
histogram(SST_prime_total_1vec,50,'normalization','pdf','displayname','6-hourly $T_o''$')
xlabel('$ T_o''$ [K]','interpreter','latex');
ylabel('Frequency','interpreter','latex');
title('DJFM, all years 2003-2018','interpreter','latex');
minTo = min([SST_prime_mean(:);SST_prime_total_1vec]);
maxTo = max([SST_prime_mean(:);SST_prime_total_1vec]);
min_label = sprintf('min = %2.2f',minTo);
max_label = sprintf('max = %2.2f',maxTo);
set(gca,'fontsize',20,'xlim',[-5 5],'xtick',[-5 -2 0 2 5],...
    'ytick',0.2:0.2:1)
set(gcf,'color','w')
legend('location','best','interpreter','latex')

update_figure_paper_size()
print(sprintf('imgs/histogram_To_prime_L_%d_%s_%d_%s',...
                L/1000,filter_type,box_num,model_str),'-dpdf')
           
fprintf('max blue To'': %f\n',max(SST_prime_mean(:)))            
fprintf('min blue To'': %f\n',min(SST_prime_mean(:)))            
fprintf('max orange To'': %f\n',max(SST_prime_total_1vec))            
fprintf('min orange To'': %f\n',min(SST_prime_total_1vec))            
            
sum(abs(SST_prime_mean(:))<2)./length(SST_prime_mean(:))
sum(abs(SST_prime_total_1vec)<2)./length(SST_prime_total_1vec)

SEM = std(SST_prime_total_1vec)/sqrt(length(SST_prime_total_1vec));  % Standard Error
zs = tinv([0.025  0.975],length(SST_prime_total_1vec)-1);            % z-Score
CI = mean(SST_prime_total_1vec) + zs*SEM;                            % Confidence Intervals


%{


[~,h] = contourf(lon_box,lat_box,nanmean(SST_box,3)');
[~,h] = contourf(lon_box,lat_box,nanmean(SST_prime,3)');
[~,h] = contourf(lon_box,lat_box,SST_time_mean_prime);
[~,h] = contourf(lon_box,lat_box,SST_time_mean_prime-nanmean(SST_prime,3)');


%}

%{

% th = text(0.01,0.05,'d)','interpreter','latex','fontsize',20,'units','normalized','BackgroundColor','w');
ax2 = subplot(3,2,1:2);
% x_plt = [year_vec fliplr(year_vec)];
% y_plt = [mean(SST_prime_total)+std(SST_prime_total) fliplr(mean(SST_prime_total)-std(SST_prime_total))]; 
% fh = fill(x_plt,y_plt,[1 1 1]*0.7); 
plot(year_vec,mean(SST_prime_total),'ks','linewidth',2,'markerfacecolor','k')
hold on
ebh1 = errorbar(year_vec,mean(SST_prime_total),std(SST_prime_total));
ebh2 = errorbar(year_vec,mean(SST_prime_total),2*std(SST_prime_total));
set(ebh1,'color','k','linewidth',2,'linestyle','none','displayname','$\sigma$')
set(ebh2,'color','k','linestyle','none','displayname','$2\sigma$')
legend([ebh1 ebh2],'location','eastoutside','interpreter','latex')
set(gca,'fontsize',20,'xlim',[2002 2019],'xtick',2003:3:2018)
xlabel('year','interpreter','latex')
ylh = title('$T_o''$ [K], DJFM','interpreter','latex');
set(gcf,'color','w','position',[1     1   720   804])
% th = text(0.01,0.1,'c)','interpreter','latex','fontsize',20,'units','normalized','BackgroundColor','w');


set(ax1,'position',[0.1300    0.0800  0.7    0.5154])
set(ax2,'position',[0.1300    0.71    0.7    0.2157])
%}