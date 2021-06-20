

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
        plot(year_vec,abCD_vec(1,:)*T_norm,'-o','color',c(1),'linewidth',2,...
            'displayname','$\alpha$');
        hold on
        yyaxis right
        plot(year_vec,abCD_vec(2,:),'color',c(2),'linewidth',2,...
            'displayname','$C_D$')
%         lh = legend('-dynamiclegend','interpreter','latex');
    case 'beta'
        h(2) = subplot(3,1,2);
        plot(year_vec,abCD_vec(1,:)*T_norm./U_norm,'--','color',c(1),'linewidth',2,...
            'displayname','$\beta$');
        hold on
        yyaxis right
        plot(year_vec,abCD_vec(2,:),'color',c(2),'linewidth',2,...
            'displayname','$C_D$')
%         lh = legend('-dynamiclegend','interpreter','latex');
    case 'alphabeta'
        h(3) = subplot(3,1,3);
        plot(year_vec,abCD_vec(1,:)*T_norm,'-o','color',c(1),'linewidth',2,...
            'displayname','$\alpha T^*/U^*$');
        hold on
        plot(year_vec,abCD_vec(2,:)*T_norm./U_norm,'--','color',c(1),'linewidth',2,...
            'displayname','$\beta T^*/U^*$')
        yyaxis right
        plot(year_vec,abCD_vec(3,:),'color',c(2),'linewidth',2,...
            'displayname','$C_D$')
        lh = legend('-dynamiclegend','interpreter','latex');
end

set(gca,'fontsize',25)
% yyaxis left
set(gca,'YColor',c(1))
yyaxis right
set(gca,'YColor',c(2))

if ~strcmp(model_str,'alphabeta')
    xlabel('')
    set(gca,'xticklabel','')
else
    set(gcf,'color','w','position',[8     1   720   804])
    xlabel('year')
    rearrange_figure(h,lh,'3x1_1_legend')
    update_figure_paper_size()
    print(sprintf('imgs/abCD_vs_year_L_%d_%s_%d_%s',...
        L/1000,filter_type,box_num,model_str),'-dpdf')
    
end









