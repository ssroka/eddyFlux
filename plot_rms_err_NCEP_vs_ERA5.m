ccc
year_vec_plt_red = [2003:2005 2010:2013 2014 2015 2018]; % red
year_vec_plt_blue = [2006:2009 2016:2017]; % blue
model_units = {'[K$^{-1}$]','[ms$^{-1}$K$^{-1}$]','[K$^{-1}$]','[ms$^{-1}$K$^{-1}$]'};
addpath ~/Documents/MATLAB/util/othercolor/
model_str_cell = {'alpha','beta','alphabeta'};
year_vec = 2003:2018;
year_vec_NCEP = 2003:2003;
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
data_base = '/Volumes/SydneySroka_Remy/eddyFlux/';
box_num = 3;
er_box_num = 3;
model_str = 'beta'; % 'alpha' 'beta' 'alphabeta'
mean_rms = zeros(length(year_vec),1);
con_str = '';
reanalysis_src = 'NCEP';% 'ERA5' 'NCEP'
% abCD_fac_vec = [-10 -8 -5 -4 -2 0.1 0.5 0.9 1 1.1 1.5 2 4 5 8 10];
abCD_fac_vec = [-10 -8 -5 -4 -2 -0.5 1 0.5 2 4 5 8 10];
% abCD_fac_vec = [0.01 0.1 1 10 100];
%%
abCD_fac_vec_1_id = find(abCD_fac_vec==1);
c_map = othercolor(296);
c_entries = round(linspace(1,size(c_map,1),length(year_vec)));
% clr_map = c_map(c_entries,:); %153 184 296
clr_map = sprintf('%s','r'*ones(1,length(year_vec)));
clr_map(ismember(year_vec,year_vec_plt_blue)) = 'b';
if strcmp(model_str,'alphabeta')
    disp_name_str{1} = '$\alpha$';
    disp_name_str{2} = '$\beta$';
else
    disp_name_str{1} = sprintf('$\\%s$',model_str);
end

j = 2;

abCD_fac_vec_1_id  = find(abCD_fac_vec==1);
for y = 1:length(year_vec)
    year = year_vec(y);
    NCEP_data = load(sprintf('rms_plot_data_%sfilt_%s_L_%d_box%d_%d_%s_%s_%d_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,'NCEP',year_vec_NCEP(1),year_vec_NCEP(end)),'abCD_fac_vec','abCD_mat','rms_er');
    ERA5_data = load(sprintf('rms_plot_data_%sfilt_%s_L_%d_box%d_%d_%s_%s_%d_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,'ERA5',year_vec(1),year_vec(end)),'abCD_fac_vec','abCD_mat','rms_er');
    h(1) = plot(abCD_fac_vec.*NCEP_data.abCD_mat(1,y),NCEP_data.rms_er(:,y),'--','linewidth',1,'color',clr_map(y),'displayname','NCEP');
    hold on
    h(2) = plot(abCD_fac_vec.*ERA5_data.abCD_mat(1,y),ERA5_data.rms_er(:,y),'-','linewidth',1,'color',clr_map(y),'displayname','ERA5');

    h(3) = plot(NCEP_data.abCD_mat(1,:),NCEP_data.rms_er(abCD_fac_vec_1_id,:),'ko','markerfacecolor','k','linewidth',2,'displayname',sprintf('RMSE at optimal %s',disp_name_str{1}));
    plot(ERA5_data.abCD_mat(1,:),ERA5_data.rms_er(abCD_fac_vec_1_id,:),'ko','markerfacecolor','k','linewidth',2,'displayname',sprintf('RMSE at optimal %s',disp_name_str{1}));

    xlabel(sprintf('$\\gamma$ %s %s %s',disp_name_str{1},model_units{j}),'interpreter','latex')
    ylabel('RMSE  $$[$$ Wm$$^{-2}]$$','interpreter','latex')
    set(gca,'fontsize',30)
    set(gcf,'color','w')
    set(gcf,'position',[1   274   996   531])
    ylim = get(gca,'ylim');
end

legend([h(1),h(2),h(3)],'interpreter','latex','location','southeast')







