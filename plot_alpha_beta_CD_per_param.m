

abCD_vec = zeros(length(abCD0),length(year_vec));

for i = 1:length(year_vec)
    year = year_vec(i);
    
    load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%s_%d',con_str,filter_type,L/1000,box_num,model_str,year_vec(i)),'abCD');
    abCD_vec(:,i) = abCD;
end

c = 'brk';

T_norm = 1;%288;
U_norm = 1;%9;

switch model_str
    case 'alpha'
        h(1) = subplot(3,1,1);
        hl1 = plot(year_vec,abCD_vec(1,:),'-o','color',c(1),'linewidth',2,...
            'displayname','$\alpha$ model');
        hold on
        h(3) = subplot(3,1,3);
        plot(year_vec,abCD_vec(2,:),'-o','color',c(1),'linewidth',2,...
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
        hl3 = plot(year_vec,abCD_vec(1,:),'-','color',c(3),'linewidth',2,...
            'displayname','$\alpha \beta$ model');
        title('$\alpha$ $[K^{-1}]$','interpreter','latex')
        hold on
        h(2) = subplot(3,1,2);
        plot(year_vec,abCD_vec(2,:),'-','color',c(3),'linewidth',2,...
            'displayname','$\beta [ms^{-1}K^{-1}]$')
        title('$\beta$ [ms$^{-1}$K$^{-1}$]','interpreter','latex')
        hold on
        h(3) = subplot(3,1,3);
        plot(year_vec,abCD_vec(3,:),'color',c(3),'linewidth',2,...
            'displayname','$C_D$')
        title('$C_D$','interpreter','latex')
        lh = legend([hl1 hl2 hl3],'interpreter','latex');
end
for i = 1:3
    h(i) = subplot(3,1,i);
    set(gca,'fontsize',25)

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








