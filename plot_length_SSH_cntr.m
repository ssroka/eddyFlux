ccc
year_vec_plt_red = [2003:2005 2010:2013 2014 2015 2018]; % red
year_vec_plt_blue = [2006:2009 2016:2017]; % blue
L        = 250000; % m
box_num = 3;
year_bounds = [2003 2018];


shift_text_vert = [zeros(1,16)];
shift_text_horiz = ones(1,16)*50;
% 2004
shift_text_vert(2) = -1.5;
shift_text_horiz(2) = -250;
% 2005
shift_text_vert(3) = 1.5;
shift_text_horiz(3) = -235;
% 2006
shift_text_vert(4) = -1.5;
% shift_text_horiz(4) = -180;
% 2008
shift_text_vert(6) = 0;
shift_text_horiz(6) = -250;
% 2010
shift_text_horiz(8) = -250;
% 2013
shift_text_vert(11) = 2.1;
shift_text_horiz(11) = -40;
% 2014
shift_text_vert(12) = -2.2;
shift_text_horiz(12) = -125;
% 2016
shift_text_horiz(14) = -250;
%{
% cntr_id = 1
shift_text_vert = [0 +0.3 -0.3 zeros(1,13)];
shift_text_horiz = ones(1,16)*50;
shift_text_horiz(10) = -200;
shift_text_horiz(12) = -200;
shift_text_vert(12) = 1;




[2003:2018;1:16]
%}
%%
cntr1 = load(sprintf('small_SSH_length_box%d_%d_%d',box_num,year_bounds(1),year_bounds(end)));
cntr2 = load(sprintf('SSH_length_box%d_%d_%d',box_num,year_bounds(1),year_bounds(end)));

cntr = [cntr1.cntr cntr2.cntr];
d = [cntr1.d;cntr2.d];
HF = cntr1.HF; % doesn't matter which
year_vec = cntr1.year_vec; % doesn't matter which
model_str = cntr1.model_str; % doesn't matter which
m = zeros(length(cntr),3); % cols: red blue all
p_corr= zeros(length(cntr),3); % cols: red blue all
d_km = d*111;
for cntr_id = 5
    close all
    %     figure(10)
    %     for i = 1:length(year_vec)
    %         plot(cntr, d(:,i),'o-','linewidth',2,'displayname',num2str(year_vec(i)))
    %         hold on
    %     end
    %     legend('-dynamiclegend')
    %     xlabel('contour height')
    %     ylabel('distance [deg]')
    %     set(gca,'fontsize',25,'xtick',cntr)
    %     set(gcf,'color','w','position',[1 215 1439 587])
    %     update_figure_paper_size()
    %     print(sprintf('imgs/ssh_length_box%d_%d_%d_%d',box_num,year_vec(1),year_vec(end),cntr_id),'-dpdf')
    
    figure(2)
    clr = 'k';
    c = 3;
    [B,BINT,R,RINT,STATS] = regress( HF(:),[ones(length(HF(:)),1) d_km(cntr_id,:)']);
    stat_vec(cntr_id,c) = STATS(3);
    % fit line
    p = polyfit(d_km(cntr_id,:),HF(:)',1);
    h_total = plot(d_km(cntr_id,:),polyval(p,d_km(cntr_id,:)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f [Wm$^{-2}$km$^{-}$], p=%2.2f',p(1),stat_vec(cntr_id,c)));
    hold on
    plot(d_km(cntr_id,:),HF(:)','x','linewidth',2,'markersize',10,'color',clr)
    for i = 1:length(year_vec)
        th = text(d_km(cntr_id,i)+0.25,HF(i),num2str(year_vec(i)));
        set(th,'fontsize',15,'color',clr,'BackgroundColor','w')
    end
    
    xlabel(sprintf('Length of %2.1fm Contour [km]',cntr(cntr_id)),'interpreter','latex')
    title('$\langle Q \rangle$ [Wm$^{-2}$]','interpreter','latex')
    set(gca,'fontsize',25)
    set(gcf,'color','w','position',[1         215        1041         587])
    lh = legend(h_total);
    set(lh,'interpreter','latex');
    m(cntr_id,c) = [1 0]*(corrcoef(d_km(cntr_id,:),HF(:)')*[0;1]);
    
    update_figure_paper_size()
    print(sprintf('imgs/ssh_length_vs_HF_box%d_%s_%d_%d_%d_allyears',box_num,model_str,year_vec(1),year_vec(end),cntr_id),'-dpdf')
    
    figure(1)
    % include p value in legend
    h_total_rb = plot(d_km(cntr_id,:),polyval(p,d_km(cntr_id,:)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f$\\times10^2$ [Wm$^{-2}$km$^{-1}$], p=%2.2f',p(1)*10^2,stat_vec(cntr_id,c)));
    % don't include p value in legend
%     h_total_rb = plot(d_km(cntr_id,:),polyval(p,d_km(cntr_id,:)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f$\\times10^2$ [Wm$^{-2}$km$^{-1}$]',p(1)*10^2));
    hold on
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
        [B,BINT,R,RINT,STATS] = regress( HF(year_ids)',[ones(length(HF(year_ids)),1) d_km(cntr_id,year_ids)']);
        stat_vec(cntr_id,c) = STATS(3);
        % fit line
        p = polyfit(d_km(cntr_id,year_ids),HF(year_ids)',1);
        % include p value in legend
        h(c) = plot(d_km(cntr_id,year_ids),polyval(p,d_km(cntr_id,year_ids)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f$\\times10^2$ [Wm$^{-2}$km$^{-1}$], p=%2.2f',p(1)*10^2,stat_vec(cntr_id,c)));
        % don't include p value in legend
%         h(c) = plot(d_km(cntr_id,year_ids),polyval(p,d_km(cntr_id,year_ids)),'linewidth',2,'color',clr,'displayname',sprintf('m = %3.2f$\\times10^2$ [Wm$^{-2}$km$^{-1}$]',p(1)*10^2));
        hold on
        plot(d_km(cntr_id,year_ids),HF(year_ids)','x','linewidth',2,'markersize',10,'color',clr)
        
        xlabel(sprintf('Length of %2.1fm Contour [km]',cntr(cntr_id)),'interpreter','latex')
        title('$\langle Q \rangle$ [Wm$^{-2}$]','interpreter','latex')
        set(gca,'fontsize',25)
        set(gcf,'color','w','position',[1         215        1041         587])
        m(cntr_id,c)= [1 0]*(corrcoef(d_km(cntr_id,year_ids),HF(year_ids)')*[0;1]);
    end
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
        for i = 1:length(year_vec_plt)
            th = text(d_km(cntr_id,year_ids(i))+shift_text_horiz(year_ids(i)),HF(year_ids(i))+shift_text_vert(year_ids(i)),num2str(year_vec_plt(i)));
            set(th,'fontsize',15,'color',clr,'BackgroundColor','w')
        end
    end

    if cntr_id == 5
            plot(8750,236,'rx','linewidth',2,'markersize',20);
            plot(8750,233,'bx','linewidth',2,'markersize',20);
            ts = text(8850,236,'stable','fontsize',25,'interpreter','latex');
            tu = text(8850,233,'unstable','fontsize',25,'interpreter','latex');
            x = [8700 9400];
            y = [231 238];
            plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'k-')
    end
    lh = legend([h h_total_rb]);
    set(lh,'interpreter','latex','location','southwest','fontsize',20);
    

    
    update_figure_paper_size()
    print(sprintf('imgs/ssh_length_vs_HF_box%d_%s_%d_%d_%d_redblue',box_num,model_str,year_vec(1),year_vec(end),cntr_id),'-dpdf')
    
    
    
    % add all year scatter to first figure
    
end

figure(3)
g(1) = plot(cntr,m(:,1),'o-','linewidth',2,'color','r','displayname','stable');
hold on
g(2) = plot(cntr,m(:,2),'o-','linewidth',2,'color','b','displayname','unstable');
g(3) = plot(cntr,m(:,3),'o-','linewidth',2,'color','k','displayname','all years');
plot(cntr,0*cntr,'--','color',[1 1 1]*0.75,'linewidth',2);
xlabel('Sea Suface Height Anomaly Contour [m]','interpreter','latex')
ylabel('correlation coefficient','interpreter','latex')
set(gca,'fontsize',25,'xtick',cntr)
set(gcf,'color','w','position',[1         215        1041         587])
lh = legend(g);
set(lh,'interpreter','latex','location','southeast','fontsize',20);

update_figure_paper_size()
print(sprintf('imgs/corr_ssh_length_vs_HF_box%d_%s_%d_%d_%d_allyears',box_num,model_str,year_vec(1),year_vec(end),cntr_id),'-dpdf')



