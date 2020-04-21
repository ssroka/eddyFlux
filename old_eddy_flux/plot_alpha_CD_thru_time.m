clear;clc;close all



load X_FFINAL_2007.mat
titles = {'$$\alpha_s^{const}$$','$$\alpha_s^{linear}$$','$$\alpha_L^{const}$$',...
          '$$\alpha_L^{linear}$$','$$C_D^{s}$$','$$C_D^{L}$$'};
      
c = 'bkgmr';
s = '*vsodp';

% titles = {'$$\alpha_s^{const}$$','$$\alpha_L^{const}$$',...
%           '$$C_D^{s}$$','$$C_D^{L}$$'};

for year = 2003:2007
%     load(sprintf('X_FFINAL_%d_const_only',year),'X');
    load(sprintf('X_FFINAL_%d',year),'X');
    aCd0 = X{year-2002};
    for i = 1:length(aCd0)
        subplot(3,2,i)
        plot(year,aCd0(i),'color','k','marker',s(i))
        hold on
        title(titles{i},'interpreter','latex')
        set(gca,'fontsize',15)
    end
end

set(gcf,'color','w','position',[277   368   725   415])
img_loc = '/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/';
print([img_loc 'model_const_thru_time'],'-dpng')
savefig([img_loc 'model_const_thru_time'])

% print([img_loc 'model_const_thru_time_const_only'],'-dpng')
% savefig([img_loc 'model_const_thru_time_const_only'])