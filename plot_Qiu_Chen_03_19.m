% clear;clc;close all

% year_vec = [2003:2018];

% data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';
% 
% addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
% addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
% addpath('/Users/ssroka/Documents/MATLAB/util/')
% addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

%% filter set up
% load('env_const.mat')
% 
% load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
%     'SST_patch','lat','lon','patch_lat','patch_lon',...
%     'slhf_patch','sshf_patch');
% L = 250000;
% box_num = 3;
% filter_type = 'fft';
% model_str = 'alpha'; % shouldn't matter, not a filter prop
for i = 1:length(year_vec)
    year = year_vec(i);
    figure(1)
    subplot(4,4,i)
    plot_HF
    hold on
    plot_SSH_contour(data_src,year,box_opt)
    title(sprintf('%d',year))
    xlabel('deg')
    ylabel('deg')
    set(gca,'fontsize',20)
    figure(2)
    subplot(4,4,i)
    plot_SST_prime(L,filter_type,model_str,box_num,lon_box,lat_box,data_src,year);
    hold on
    plot_SSH_contour(data_src,year,box_opt)
    title(sprintf('%d',year))
    xlabel('deg')
    ylabel('deg')
    set(gca,'fontsize',20)
end

for j = 1:2
figure(j)
set(gcf,'position',[45          38        1354         767],'color','w')
update_figure_paper_size()
if j==1
print(sprintf('imgs/Qiu_Chen_repeat_HF_%d_%d',year_vec(1),year_vec(end)),'-dpdf')
% print(sprintf('imgs/HF_and_SSH_%d_%d',year_vec(1),year_vec(end)),'-dpdf')
else 
    print(sprintf('imgs/Qiu_Chen_repeat_SST_Prime_%d_%d',year_vec(1),year_vec(end)),'-dpdf')
end
% print(sprintf('imgs/HF_and_SSH_%d_%d',year_vec(1),year_vec(end)),'-dpdf')
end