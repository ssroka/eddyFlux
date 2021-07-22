ccc
year_vec_plt_red = [2003:2005 2010:2013 2014 2015 2018]; % red
year_vec_plt_blue = [2006:2009 2016:2017]; % blue
L        = 250000; % m
box_num = 3;
year_bounds = [2003 2018];
%%
load(sprintf('SSH_length_box%d_%d_%d',box_num,year_bounds(1),year_bounds(end)))
m = zeros(length(cntr),3); % cols: red blue all
for cntr_id = 1:length(cntr)
    close all
    figure(10)
    for i = 1:length(year_vec)
        plot(cntr, d(:,i),'o-','linewidth',2,'displayname',num2str(year_vec(i)))
        hold on
    end
    legend('-dynamiclegend')
    xlabel('contour height')
    ylabel('distance [deg]')
    set(gca,'fontsize',25,'xtick',cntr)
    set(gcf,'color','w','position',[1 215 1439 587])
    update_figure_paper_size()
    print(sprintf('imgs/ssh_length_box%d_%d_%d_%d',box_num,year_vec(1),year_vec(end),cntr_id),'-dpdf')
    for c = 1:2
        if c == 1
            year_vec_plt = year_vec_plt_red;
            year_ids = find(ismember(year_bounds(1):year_bounds(2),year_vec_plt))';
            clr = 'r';
        else
            year_vec_plt = year_vec_plt_blue;
            year_ids = find(ismember(year_bounds(1):year_bounds(2),year_vec_plt))';
            clr = 'b';
        end
        figure(1)
        [B,BINT,R,RINT,STATS] = regress( HF(year_ids)',[ones(length(HF(year_ids)),1) d(cntr_id,year_ids)']);
        stat_vec(cntr_id,c) = STATS(3);
        plot(d(cntr_id,year_ids),HF(year_ids)','x','linewidth',2,'markersize',10,'color',clr)
        hold on
        for i = 1:length(year_vec_plt)
            th = text(d(cntr_id,year_ids(i))+0.25,HF(year_ids(i)),num2str(year_vec_plt(i)));
            set(th,'fontsize',15,'color',clr)
        end
        
        % fit line
        p = polyfit(d(cntr_id,year_ids),HF(year_ids)',1);
        h(c) = plot(d(cntr_id,year_ids),polyval(p,d(cntr_id,year_ids)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f $W/m^2/deg$, p=%2.2f',p(1),stat_vec(cntr_id,c)));
        xlabel(sprintf('SSH length of %2.1fm contour [deg]',cntr(cntr_id)),'interpreter','latex')
        ylabel('mean total flux $[W/m^2]$','interpreter','latex')
        set(gca,'fontsize',25)
        set(gcf,'color','w','position',[1 215 1439 587])
        m(cntr_id,c) = [1 0]*(corrcoef(d(cntr_id,year_ids),HF(year_ids)')*[0;1]);
    end
    lh = legend(h);
    set(lh,'interpreter','latex');
    
    update_figure_paper_size()
    print(sprintf('imgs/ssh_length_vs_HF_box%d_%s_%d_%d_%d_redblue',box_num,model_str,year_vec(1),year_vec(end),cntr_id),'-dpdf')
    
    
    figure(2)
    clr = 'k';
    c = 1;
    [B,BINT,R,RINT,STATS] = regress( HF(:),[ones(length(HF(:)),1) d(cntr_id,:)']);
    stat_vec(cntr_id,c) = STATS(3);
    plot(d(cntr_id,:),HF(:)','x','linewidth',2,'markersize',10,'color',clr)
    hold on
    for i = 1:length(year_vec)
        th = text(d(cntr_id,i)+0.25,HF(i),num2str(year_vec(i)));
        set(th,'fontsize',15,'color',clr)
    end
    
    % fit line
    p = polyfit(d(cntr_id,:),HF(:)',1);
    h_total = plot(d(cntr_id,:),polyval(p,d(cntr_id,:)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f $W/m^2/deg$, p=%2.2f',p(1),stat_vec(cntr_id,c)));
    xlabel(sprintf('SSH length of %2.1fm contour [deg]',cntr(cntr_id)),'interpreter','latex')
    ylabel('mean total flux $[W/m^2]$','interpreter','latex')
    set(gca,'fontsize',25)
    set(gcf,'color','w','position',[1 215 1439 587])
    lh = legend(h_total);
    set(lh,'interpreter','latex');
    m(cntr_id,3) = [1 0]*(corrcoef(d(cntr_id,:),HF(:)')*[0;1]);
    
    update_figure_paper_size()
    print(sprintf('imgs/ssh_length_vs_HF_box%d_%s_%d_%d_%d_allyears',box_num,model_str,year_vec(1),year_vec(end),cntr_id),'-dpdf')
    
end

close all
g(1) = plot(cntr,m(:,1),'o-','linewidth',2,'color','r','displayname','stable');
hold on
g(2) = plot(cntr,m(:,2),'o-','linewidth',2,'color','b','displayname','unstable');
g(3) = plot(cntr,m(:,3),'o-','linewidth',2,'color','k','displayname','all years');
plot(cntr,0*cntr,'--','color',[1 1 1]*0.75,'linewidth',2);
xlabel('Sea Suface Height Anomaly Contour [m]','interpreter','latex')
ylabel('correlation coefficient','interpreter','latex')
set(gca,'fontsize',25,'xtick',cntr)
set(gcf,'color','w','position',[1 215 1439 587])
lh = legend(g);
set(lh,'interpreter','latex','location','southeast');

update_figure_paper_size()
print(sprintf('imgs/corr_ssh_length_vs_HF_box%d_%s_%d_%d_%d_allyears',box_num,model_str,year_vec(1),year_vec(end),cntr_id),'-dpdf')



