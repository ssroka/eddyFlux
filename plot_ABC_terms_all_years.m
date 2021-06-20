
if strcmp(model_str,'beta')
    num_comp = 4;
else
    num_comp = 8;
end

mean_comp = zeros(num_comp,length(year_vec));

for i = 1:length(year_vec)
    year = year_vec(i);
    load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
    switch model_str
        case 'alpha'
            mean_comp(1,i) = sum(sum(nanmean(A,3)));
            mean_comp(2,i) = sum(sum(nanmean(B,3)));
            mean_comp(3,i) = sum(sum(nanmean(C1,3)));
            mean_comp(4,i) = sum(sum(nanmean(C2,3)));
            mean_comp(5,i) = sum(sum(nanmean(C3,3)));
            mean_comp(6,i) = sum(sum(nanmean(D,3)));
            mean_comp(7,i) = sum(sum(nanmean(E1,3)));
            mean_comp(8,i) = sum(sum(nanmean(E2,3)));
            comp_str = {'A','B','C1','C2','C3','D','E1','E2'};
            
        case 'beta'
            mean_comp(1,i) = sum(sum(nanmean(Ab,3)));
            mean_comp(2,i) = sum(sum(nanmean(Bb,3)));
            mean_comp(3,i) = sum(sum(nanmean(Cb,3)));
            mean_comp(4,i) = sum(sum(nanmean(Db,3)));
            comp_str = {'A','B','C','D'};
        case 'alphabeta'
            mean_comp(1,i) = sum(sum(nanmean(Aab,3)));
            mean_comp(2,i) = sum(sum(nanmean(Bab,3)));
            mean_comp(3,i) = sum(sum(nanmean(C1ab,3)));
            mean_comp(4,i) = sum(sum(nanmean(C2ab,3)));
            mean_comp(5,i) = sum(sum(nanmean(D1ab,3)));
            mean_comp(6,i) = sum(sum(nanmean(D2ab,3)));
            mean_comp(7,i) = sum(sum(nanmean(E1ab,3)));
            mean_comp(8,i) = sum(sum(nanmean(E2ab,3)));
            comp_str = {'A','B','C1','C2','D1','D2','E1','E2'};
    end
    
end
mean_comp = mean_comp/(m*n);

c = 'brk';
switch model_str
    case 'alpha'
        [titles_cell] = term_names(model_str);
        for i = 1:num_comp
            
            if i > 1
                ax(2) = subplot(3,2,2);
            else
                ax(1) = subplot(3,2,1);
            end
            
            h(1) = plot(year_vec,mean_comp(i,:),'-','linewidth',2,'displayname',titles_cell{i});
            set(gca,'fontsize',15)
            hold on
            
            if i > 1
                set(gca,'ylim',[-0.2 0.8])
                title('$\alpha$ model, terms B-E [Wm$^{-2}$]','interpreter','latex')
                lh(1) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
            else
                
                set(gca,'ylim',[250 300])
                title(['$\alpha$ model, term ' titles_cell{1} ' [Wm$^{-2}$]'],'interpreter','latex')
            end
            xlabel('year')
        end
        set(gcf,'color','w','position',[61 221 1359 581])
        
    case 'beta'
        %         figure(2)
        [titles_cell] = term_names(model_str);
        
        for i = 1:num_comp
            if i > 1
                ax(4) = subplot(3,2,4);
            else
                ax(3) = subplot(3,2,3);
            end
            
            
            h(2) = plot(year_vec,mean_comp(i,:),'-','linewidth',2,'displayname',titles_cell{i});
            
            set(gca,'fontsize',15)
            hold on
            if i > 1
                set(gca,'ylim',[-0.2 0.8])
                lh(2) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
                title('$\beta$ model, terms B-D [Wm$^{-2}$]','interpreter','latex')
                
            else
                set(gca,'ylim',[250 300])
                title(['$\beta$ model, term ' titles_cell{1} ' [Wm$^{-2}$]'],'interpreter','latex')
                
            end
            xlabel('year')
            
        end
        
        set(gcf,'color','w','position',[61 221 1359 581])
        
    case 'alphabeta'
        [titles_cell] = term_names(model_str);
        
        for i = 1:num_comp
            
            if i > 1
                ax(6) = subplot(3,2,6);
            else
                ax(5) = subplot(3,2,5);
            end
            
            h(3) = plot(year_vec,mean_comp(i,:),'linewidth',2,'displayname',titles_cell{i});
            set(gca,'fontsize',15)
            hold on
            if i > 1
                set(gca,'ylim',[-0.2 0.8])
                title('$\alpha \beta$ model, terms B-E [Wm$^{-2}$]','interpreter','latex')
                lh(3) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
            else
                set(gca,'ylim',[250 300])
                title(['$\alpha \beta$ model, term ' titles_cell{1} ' [Wm$^{-2}$]' ],'interpreter','latex')
            end
            xlabel('year')
        end
        rearrange_figure(ax,lh,'3x2_3_legend');
        %         title('$\alpha \beta$ model','interpreter','latex')
        %         format_fig(h,ax)
        set(gcf,'color','w','position',[61 221 1359 581])
        update_figure_paper_size()
        print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_all_years_L_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
        
        
end







function [titles_cell] = term_names(model_str)

switch model_str
    case 'alpha'
        titles_cell = {'A:$\rho_a C_D\overline{\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            'B:$\rho_a C_D\overline{U''\Delta h''}$',...
            'C1:$\rho_a C_D\overline{\overline{U}T_o''\alpha \Delta h''}$',...
            'C2:$\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$',...
            'C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$',...
            'D:$\rho_a C_D\overline{U''\overline{\Delta h}}$',...
            'E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$',...
            'E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$'};
        
        
    case 'beta'
        
        titles_cell = {'A:$\overline{C_D\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            'B:$C_D\overline{\overline{U}\Delta h''}$',...
            'C:$C_D\overline{\beta T_o'' \Delta h''}$',...
            'D:$C_D\overline{\beta T_o''\overline{\Delta h}}$'};
        
        
    case 'alphabeta'
        
        titles_cell = {'A:$\rho_a C_D\overline{\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            'B:$\rho_a C_D\overline{\overline{U}\Delta h''}$',...
            'C1:$\rho_a C_D\overline{\overline{U}T_o''\alpha \overline{\Delta h}}$',...
            'C2:$\rho_a C_D\overline{\overline{U}T_o''\alpha \Delta h''}$',...
            'D1:$\rho_a C_D\overline{T_o''\beta \overline{\Delta h}}$',...
            'D2:$\rho_a C_D\overline{T_o''\beta  \Delta h''}$',...
            'E1:$\rho_a C_D\overline{ (T_o'')^2\alpha\beta \overline{\Delta h}}$',...
            'E2:$\rho_a C_D\overline{ (T_o'')^2\alpha\beta \Delta h''}$'};
        
end

end










