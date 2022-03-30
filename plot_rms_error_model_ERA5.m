ccc
year_vec_plt_red = [2003:2005 2010:2013 2014 2015 2018]; % red
year_vec_plt_blue = [2006:2009 2016:2017]; % blue
model_units = {'[K$^{-1}$]','[ms$^{-1}$K$^{-1}$]','[K$^{-1}$]','[ms$^{-1}$K$^{-1}$]'};
addpath ~/Documents/MATLAB/util/othercolor/
model_str_cell = {'alpha','beta','alphabeta'};
year_vec = 2003;
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
data_base = '/Volumes/SydneySroka_Remy/eddyFlux/';
box_num = 3;
er_box_num = 3;
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
for j = 2
    rms_er = zeros(length(abCD_fac_vec),length(year_vec));
    model_str = model_str_cell{j};
    for i = 1:length(abCD_fac_vec)
        if strcmp(model_str,'alphabeta')
            abCD_factor = [abCD_fac_vec(i) abCD_fac_vec(i) 1];
            abCD_mat = zeros(3,length(year_vec));
            disp_name_str{1} = '$\alpha$';
            disp_name_str{2} = '$\beta$';
            R2 = zeros(length(year_vec),2);
            
        else
            abCD_factor = [abCD_fac_vec(i) 1];
            abCD_mat = zeros(2,length(year_vec));
            disp_name_str{1} = sprintf('$\\%s$',model_str);
            R2 = zeros(length(year_vec),1);
            
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
        for y = 1:length(year_vec)
            year = year_vec(y);
            load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%d_%s_%s_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,reanalysis_src,year),'abCD');
            load(sprintf('mean_rms_err_%s_%s_%d%s',model_str,reanalysis_src,year,abCD_fac_str),'mean_rms')
            rms_er(i,y) = mean_rms;
            abCD_mat(:,y) = abCD(:);
        end
        
    end
    
    save(sprintf('rms_plot_data_%sfilt_%s_L_%d_box%d_%d_%s_%s_%d_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,reanalysis_src,year_vec(1),year_vec(end)),'abCD_fac_vec','abCD_mat','rms_er')
    
    figure(10+j)
    for y = 1:length(year_vec)
        h(1) = plot(abCD_fac_vec.*abCD_mat(1,y),rms_er(:,y),'o-','linewidth',1,'color',clr_map(y),'displayname',disp_name_str{1});
        [B,BINT,R,RINT,STATS] = regress(rms_er(:,y),[ones(length(abCD_fac_vec.*abCD_mat(1,y)),1) abCD_fac_vec'.*abCD_mat(1,y) (abCD_fac_vec'.*abCD_mat(1,y)).^2]);
        R2(y,1) = STATS(1);
        mid(y,1) = -B(2)/(2*B(3));
        hold on
    end
    h(2) = plot(abCD_mat(1,:),rms_er(abCD_fac_vec_1_id,:),'ko','markerfacecolor','k','linewidth',2,'displayname',sprintf('RMSE at optimal %s',disp_name_str{1}));
    legend([h(2)],'interpreter','latex')
    xlabel(sprintf('$\\gamma$ %s %s %s',disp_name_str{1},model_units{j}),'interpreter','latex')
    ylabel('RMSE  $$[$$ Wm$$^{-2}]$$','interpreter','latex')
    set(gca,'fontsize',30)
    set(gcf,'color','w')
    set(gcf,'position',[1   274   996   531])
    ylim = get(gca,'ylim');
    drawnow
    pause(1)% don't remove not saving properly without these
    update_figure_paper_size()
    pause(1)% don't remove not saving properly without these
    print(sprintf('%simgs/rms_error_%d_%s_box%d_%s_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,reanalysis_src,year,abCD_fac_str),'-dpdf')
    
    if strcmp(model_str,'alphabeta')
        figure(20+j)
        for y = 1:length(year_vec)
            h(3) = plot(abCD_fac_vec.*abCD_mat(2,y),rms_er(:,y),'-o','linewidth',1,'color',clr_map(y),'displayname',disp_name_str{2});
            [B,BINT,R,RINT,STATS] = regress(rms_er(:,y),[ones(length(abCD_fac_vec.*abCD_mat(2,y)),1) abCD_fac_vec'.*abCD_mat(2,y) (abCD_fac_vec'.*abCD_mat(2,y)).^2]);
            R2(y,2) = STATS(1);
            mid(y,2) = -B(2)/(2*B(3));
            hold on
        end
        h(4) = plot(abCD_mat(2,:),rms_er(abCD_fac_vec_1_id,:),'ko','markerfacecolor','k','linewidth',2,'displayname',sprintf('RMSE at optimal %s',disp_name_str{2}));
        legend([h(4)],'interpreter','latex')
        xlabel(sprintf('$\\gamma$ %s %s %s',disp_name_str{2},model_units{j+1}),'interpreter','latex')
        ylabel('RMSE  $$[$$ Wm$$^{-2}]$$','interpreter','latex')
        set(gca,'fontsize',30)
        set(gcf,'color','w')
        set(gcf,'position',[1   274   996   531])
        ylim = get(gca,'ylim');
        drawnow
    pause(1)% don't remove not saving properly without these
    update_figure_paper_size()
    pause(1)% don't remove not saving properly without these
    print(sprintf('%simgs/rms_error_%d_%s_box%d_%s_beta_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,reanalysis_src,year,abCD_fac_str),'-dpdf')
    
    end
  
    
    figure(j)
    for i = 1:length(abCD_fac_vec)
        mean_er_all_yrs = mean(rms_er(i,:));
        mean_param_all_yrs = mean(abCD_fac_vec(i).*abCD_mat(1,:));
        er_bar_x = std(abCD_fac_vec(i).*abCD_mat(1,:));
        er_bar_y = std(rms_er(i,:));
        h(1) = plot(mean_param_all_yrs,mean_er_all_yrs,'ks','markersize',10,'linewidth',2,'displayname',disp_name_str{1});
        hold on
        errorbar(mean_param_all_yrs,mean_er_all_yrs,er_bar_y,er_bar_y,er_bar_x,er_bar_x,'k','linewidth',2)
    end
    legend([h(1)],'interpreter','latex')
    if strcmp(model_str,'alphabeta')
        for i = 1:length(abCD_fac_vec)
            mean_er_all_yrs = mean(rms_er(i,:));
            mean_param_all_yrs = mean(abCD_fac_vec(i).*abCD_mat(2,:));
            er_bar_x = std(abCD_fac_vec(i).*abCD_mat(2,:));
            er_bar_y = std(rms_er(i,:));
            h(2) = plot(mean_param_all_yrs,mean_er_all_yrs,'bo','markersize',10,'linewidth',2,'displayname',disp_name_str{2});
            hold on
            errorbar(mean_param_all_yrs,mean_er_all_yrs,er_bar_y,er_bar_y,er_bar_x,er_bar_x,'b','linewidth',2)
        end
        legend([h(1:2)],'interpreter','latex')
    end
    xlabel(model_str,'interpreter','latex')
    ylabel('rms error')
    set(gca,'fontsize',20,'ylim',ylim)
    set(gcf,'color','w')
    set(gcf,'position',[1   274   996   531])
    title('the mean and standard deviation across all years')
    update_figure_paper_size()
    print(sprintf('%simgs/rms_error_erbar_%d_%s_box%d_%s_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,reanalysis_src,year,abCD_fac_str),'-dpdf')
    
end

%     legend('\alpha','\beta','\alpha\beta','location','northwest')







