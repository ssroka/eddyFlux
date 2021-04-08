clear;clc;close all

year_vec = [2003:2018];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

%% filter set up
load('env_const.mat')

load(sprintf('%sERA5_patch_data_%d.mat',data_src,2003),...
    'SST_patch','lat','lon','patch_lat','patch_lon',...
    'slhf_patch','sshf_patch');

for i = 1:length(year_vec)
    year = year_vec(i);
    subplot(4,4,i)
    plot_HF
    hold on
    plot_SSH_contour
    title(sprintf('%d',year))
    xlabel('deg')
    ylabel('deg')
    set(gca,'fontsize',20)
end


set(gcf,'position',[45          38        1354         767],'color','w')

update_figure_paper_size()
print(sprintf('imgs/Qiu_Chen_repeat_%d_%d',year_vec(1),year_vec(end)),'-dpdf')
print(sprintf('imgs/HF_and_SSH_%d_%d',year_vec(1),year_vec(end)),'-dpdf')