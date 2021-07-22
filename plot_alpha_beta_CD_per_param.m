<<<<<<< HEAD
% year_vec_plt_red = [2003:2005 2010:2013 2014 2015 2018]; % red
% year_vec_plt_blue = [2006:2009 2016:2017]; % blue

year_vec_plt_red = [2003 2005.5;
    2009.5 2015.5;
    2017.5 2018.5]; % red
year_vec_plt_blue = [2005.5 2009.5;
    2015.5 2017.5]; % blue
=======

>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28

abCD_vec = zeros(length(abCD0),length(year_vec));

for i = 1:length(year_vec)
    year = year_vec(i);
<<<<<<< HEAD
=======
    
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
    load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%s_%d',con_str,filter_type,L/1000,box_num,model_str,year_vec(i)),'abCD');
    abCD_vec(:,i) = abCD;
end

<<<<<<< HEAD
c = 'kkk';
=======
c = 'brk';
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28

T_norm = 1;%288;
U_norm = 1;%9;

switch model_str
    case 'alpha'
        h(1) = subplot(3,1,1);
<<<<<<< HEAD
        hl1 = plot(year_vec,abCD_vec(1,:),'-','color',c(1),'linewidth',2,...
            'displayname','$\alpha$ model');
        hold on
        h(3) = subplot(3,1,3);
        plot(year_vec,abCD_vec(2,:),'-','color',c(1),'linewidth',2,...
=======
        hl1 = plot(year_vec,abCD_vec(1,:),'-o','color',c(1),'linewidth',2,...
            'displayname','$\alpha$ model');
        hold on
        h(3) = subplot(3,1,3);
        plot(year_vec,abCD_vec(2,:),'-o','color',c(1),'linewidth',2,...
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            'displayname','$C_D$')
        hold on
        %         lh = legend('-dynamiclegend','interpreter','latex');
    case 'beta'
        h(2) = subplot(3,1,2);
        hl2 = plot(year_vec,abCD_vec(1,:),'--','color',c(2),'linewidth',2,...
            'displayname','$\beta$ model');
        hold on
        h(3) = subplot(3,1,3);
        plot(year_vec,abCD_vec(2,:),'--','color',c(2),'linewidth',2,...
            'displayname','$C_D$')
        hold on
        %         lh = legend('-dynamiclegend','interpreter','latex');
    case 'alphabeta'
        h(1) = subplot(3,1,1);
<<<<<<< HEAD
        hl3 = plot(year_vec,abCD_vec(1,:),':','color',c(3),'linewidth',2,...
=======
        hl3 = plot(year_vec,abCD_vec(1,:),'-','color',c(3),'linewidth',2,...
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            'displayname','$\alpha \beta$ model');
        title('$\alpha$ $[K^{-1}]$','interpreter','latex')
        hold on
        h(2) = subplot(3,1,2);
<<<<<<< HEAD
        plot(year_vec,abCD_vec(2,:),':','color',c(3),'linewidth',2,...
=======
        plot(year_vec,abCD_vec(2,:),'-','color',c(3),'linewidth',2,...
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            'displayname','$\beta [ms^{-1}K^{-1}]$')
        title('$\beta$ [ms$^{-1}$K$^{-1}$]','interpreter','latex')
        hold on
        h(3) = subplot(3,1,3);
<<<<<<< HEAD
        plot(year_vec,abCD_vec(3,:),':','color',c(3),'linewidth',2,...
=======
        plot(year_vec,abCD_vec(3,:),'color',c(3),'linewidth',2,...
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            'displayname','$C_D$')
        title('$C_D$','interpreter','latex')
        lh = legend([hl1 hl2 hl3],'interpreter','latex');
end
for i = 1:3
    h(i) = subplot(3,1,i);
    set(gca,'fontsize',25)
<<<<<<< HEAD
    if i==1
        set(gca,'ylim',[-1 1]*0.02)
    elseif i==2
        set(gca,'ylim',[0.125 0.4])
    else
        set(gca,'ylim',[1.39 1.45]*1e-3)
    end
    
    
    if i==1 || i==2
        xlabel('','interpreter','latex')
        ylabel('','interpreter','latex')
        set(gca,'xticklabel','')
    elseif i==3 && strcmp(model_str,'alphabeta')
        h(3) = subplot(3,1,3);
        ylabel('','interpreter','latex')
        xlabel('year','interpreter','latex')
        set(gca,'xtick',unique([year_vec(1):3:year_vec(end) year_vec(end)]))
        set(gcf,'color','w','position',[8   214   870   591])
        rearrange_figure(h,lh,'3x1_1_legend')
        
    end
    if strcmp(model_str,'alphabeta') % avoid running more than once
        for j = 1:size(year_vec_plt_red,1)
            x1 = year_vec_plt_red(j,1);
            x2 = year_vec_plt_red(j,2);
            y1 = get(gca,'ylim')*[1;0];
            y2 = get(gca,'ylim')*[0;1];
            fill([x1 x2 x2 x1],[y1 y1 y2 y2],'r','facealpha',0.15,'LineStyle','none')
        end
        for j = 1:size(year_vec_plt_blue,1)
            x1 = year_vec_plt_blue(j,1);
            x2 = year_vec_plt_blue(j,2);
            y1 = get(gca,'ylim')*[1;0];
            y2 = get(gca,'ylim')*[0;1];
            fill([x1 x2 x2 x1],[y1 y1 y2 y2],'b','facealpha',0.15,'LineStyle','none')
        end
        lh = legend([hl1 hl2 hl3],'interpreter','latex');
        if i==3
            rearrange_figure(h,lh,'3x1_1_legend')
            update_figure_paper_size()
            print(sprintf('imgs/abCD_vs_year_L_%d_%s_%d_%s',...
                L/1000,filter_type,box_num,'per_param_abCD'),'-dpdf')
            
        end
    end
    
end
=======

if i==1 || i==2
    xlabel('')
    set(gca,'xticklabel','')
elseif i==3 && strcmp(model_str,'alphabeta')
    h(3) = subplot(3,1,3);
    xlabel('year')
    set(gcf,'color','w','position',[8     1   720   804])
    rearrange_figure(h,lh,'3x1_1_legend')
    update_figure_paper_size()
    print(sprintf('imgs/abCD_vs_year_L_%d_%s_%d_%s',...
        L/1000,filter_type,box_num,'per_param_abCD'),'-dpdf')
    
end
end
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28








