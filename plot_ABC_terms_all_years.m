
if strcmp(model_str,'beta')
    num_comp = 4;
else
    num_comp = 8;
end

<<<<<<< HEAD
c_map = othercolor(184);
c_entries = round(linspace(1,size(c_map,1),num_comp));
clr_map = c_map(c_entries,:); %153 184 296

=======
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
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
<<<<<<< HEAD
                ax(2) = subplot(4,1,2);
            else
                continue
            end
            
            h(1) = plot(year_vec,mean_comp(i,:),'-','color',clr_map(i,:),'linewidth',2,'displayname',titles_cell{i});
=======
                ax(2) = subplot(3,2,2);
            else
                ax(1) = subplot(3,2,1);
            end
            
            h(1) = plot(year_vec,mean_comp(i,:),'-','linewidth',2,'displayname',titles_cell{i});
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            set(gca,'fontsize',15)
            hold on
            
            if i > 1
<<<<<<< HEAD
                set(gca,'ylim',[-0.4 1])
                title('$\alpha$ model $\left(\overline{Q^\alpha_2}\right)_{\overline{t}}$ through $\left(\overline{Q^\alpha_4}\right)_{\overline{t}}$ [Wm$^{-2}$]','interpreter','latex')
                lh(2) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
            else
                
                %                 set(gca,'ylim',[250 300])
                %                 title(['$\alpha$ model $\left(\overline{Q^\alpha_1}\right)_{\overline{t}}$ [Wm$^{-2}$]'],'interpreter','latex')
                % %                 th = text(1,1.3,'$\alpha$ model','units','normalized','fontsize',20,'interpreter','latex');
                
            end
            xlabel('year')
        end
        set(gcf,'color','w','position',[1     1   720   804])
=======
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
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
        
    case 'beta'
        %         figure(2)
        [titles_cell] = term_names(model_str);
        
        for i = 1:num_comp
            if i > 1
<<<<<<< HEAD
                ax(3) = subplot(4,1,3);
            else
                continue
            end
            
            
            h(2) = plot(year_vec,mean_comp(i,:),'-','color',clr_map(i,:),'linewidth',2,'displayname',titles_cell{i});
=======
                ax(4) = subplot(3,2,4);
            else
                ax(3) = subplot(3,2,3);
            end
            
            
            h(2) = plot(year_vec,mean_comp(i,:),'-','linewidth',2,'displayname',titles_cell{i});
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            
            set(gca,'fontsize',15)
            hold on
            if i > 1
<<<<<<< HEAD
                set(gca,'ylim',[-0.4 1])
                lh(3) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
                title('$\beta$ model $\left(\overline{Q^\beta_2}\right)_{\overline{t}}$ through $\left(\overline{Q^\beta_4}\right)_{\overline{t}}$ [Wm$^{-2}$]','interpreter','latex')
                
            else
                %                 set(gca,'ylim',[250 300])
                %                 title(['$\left(\overline{Q^\beta_1}\right)_{\overline{t}}$ [Wm$^{-2}$]'],'interpreter','latex')
                %                 th = text(1,1.3,'$\beta$ model','units','normalized','fontsize',20,'interpreter','latex');
=======
                set(gca,'ylim',[-0.2 0.8])
                lh(2) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
                title('$\beta$ model, terms B-D [Wm$^{-2}$]','interpreter','latex')
                
            else
                set(gca,'ylim',[250 300])
                title(['$\beta$ model, term ' titles_cell{1} ' [Wm$^{-2}$]'],'interpreter','latex')
                
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
            end
            xlabel('year')
            
        end
        
<<<<<<< HEAD
        set(gcf,'color','w','position',[1     1   720   804])
=======
        set(gcf,'color','w','position',[61 221 1359 581])
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
        
    case 'alphabeta'
        [titles_cell] = term_names(model_str);
        
        for i = 1:num_comp
            
            if i > 1
<<<<<<< HEAD
                ax(4) = subplot(4,1,4);
            else
                ax(1) = subplot(4,1,1);
            end
            
            h(3) = plot(year_vec,mean_comp(i,:),'color',clr_map(i,:),'linewidth',2,'displayname',titles_cell{i});
            set(gca,'fontsize',15)
            hold on
            if i > 1
                set(gca,'ylim',[-0.4 1])
                title('$\alpha\beta$ model $\left(\overline{Q^{\alpha\beta}_2}\right)_{\overline{t}}$ through $\left(\overline{Q^{\alpha\beta}_8}\right)_{\overline{t}}$ [Wm$^{-2}$]','interpreter','latex')
                lh(4) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
            else
                set(gca,'ylim',[250 300])
                title('$\left(\overline{Q^\alpha_1}\right)_{\overline{t}}=\left(\overline{Q^\beta_1}\right)_{\overline{t}}=\left(\overline{Q^{\alpha\beta}_1}\right)_{\overline{t}}$ [Wm$^{-2}$]','interpreter','latex')
                %                 th = text(1,1.3,'$\alpha\beta$ model','units','normalized','fontsize',20,'interpreter','latex');
                lh(1) = legend('-dynamiclegend','interpreter','latex','location','eastoutside');
            end
            xlabel('year')
        end
        set(gcf,'color','w','position',[1     1   720   804])
        
        rearrange_figure(ax,lh,'4x1_4_legend');
        %         title('$\alpha \beta$ model','interpreter','latex')
        %         format_fig(h,ax)
        update_figure_paper_size()
        print(sprintf('%simgs/ABC_all_years_L_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
=======
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
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
        
        
end







function [titles_cell] = term_names(model_str)

switch model_str
    case 'alpha'
<<<<<<< HEAD
        titles_cell = {...
            '$\left(\overline{Q^\alpha_1}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            '$\left(\overline{Q^\alpha_2}\right)_{\overline{t}}$: $\overline{\rho_a C_DU''\Delta h''}$',...
            '$\left(\overline{Q^\alpha_3}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}T_o''\alpha \Delta h''}$',...
            '$\left(\overline{Q^\alpha_4}\right)_{\overline{t}}$: $\overline{\rho_a C_DU''T_o''\overline{\alpha \Delta h}}$',...
            '$\left(\overline{Q^\alpha_5}\right)_{\overline{t}}$: $\overline{\rho_a C_DU'' T_o''\alpha \Delta h''}$',...
            '$\left(\overline{Q^\alpha_6}\right)_{\overline{t}}$: $\overline{\rho_a C_DU''\overline{\Delta h}}$',...
            '$\left(\overline{Q^\alpha_7}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U} T_o''\overline{\alpha \Delta h}}$',...
            '$\left(\overline{Q^\alpha_8}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}  \Delta h''}$'};
=======
        titles_cell = {'A:$\rho_a C_D\overline{\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            'B:$\rho_a C_D\overline{U''\Delta h''}$',...
            'C1:$\rho_a C_D\overline{\overline{U}T_o''\alpha \Delta h''}$',...
            'C2:$\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$',...
            'C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$',...
            'D:$\rho_a C_D\overline{U''\overline{\Delta h}}$',...
            'E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$',...
            'E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$'};
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
        
        
    case 'beta'
        
<<<<<<< HEAD
        titles_cell = {'$\left(\overline{Q^\beta_1}\right)_{\overline{t}}$: $\overline{\rho_aC_D\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            '$\left(\overline{Q^\beta_2}\right)_{\overline{t}}$: $\overline{\rho_aC_D\overline{U}\Delta h''}$',...
            '$\left(\overline{Q^\beta_3}\right)_{\overline{t}}$: $\overline{\rho_aC_D\beta T_o'' \Delta h''}$',...
            '$\left(\overline{Q^\beta_4}\right)_{\overline{t}}$: $\overline{\rho_aC_D\beta T_o''\overline{\Delta h}}$'};
=======
        titles_cell = {'A:$\overline{C_D\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            'B:$C_D\overline{\overline{U}\Delta h''}$',...
            'C:$C_D\overline{\beta T_o'' \Delta h''}$',...
            'D:$C_D\overline{\beta T_o''\overline{\Delta h}}$'};
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
        
        
    case 'alphabeta'
        
<<<<<<< HEAD
        titles_cell = {...
            '$\left(\overline{Q^{\alpha,\beta,\alpha\beta}_1}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            '$\left(\overline{Q^{\alpha\beta}_2}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}\Delta h''}$',...
            '$\left(\overline{Q^{\alpha\beta}_3}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}T_o''\alpha \overline{\Delta h}}$',...
            '$\left(\overline{Q^{\alpha\beta}_4}\right)_{\overline{t}}$: $\overline{\rho_a C_D\overline{U}T_o''\alpha \Delta h''}$',...
            '$\left(\overline{Q^{\alpha\beta}_5}\right)_{\overline{t}}$: $\overline{\rho_a C_DT_o''\beta \overline{\Delta h}}$',...
            '$\left(\overline{Q^{\alpha\beta}_6}\right)_{\overline{t}}$: $\overline{\rho_a C_DT_o''\beta  \Delta h''}$',...
            '$\left(\overline{Q^{\alpha\beta}_7}\right)_{\overline{t}}$: $\overline{\rho_a C_D (T_o'')^2\alpha\beta \overline{\Delta h}}$',...
            '$\left(\overline{Q^{\alpha\beta}_8}\right)_{\overline{t}}$: $\overline{\rho_a C_D (T_o'')^2\alpha\beta \Delta h''}$'};
=======
        titles_cell = {'A:$\rho_a C_D\overline{\overline{U}\hspace{1mm}\overline{\Delta h}}$',...
            'B:$\rho_a C_D\overline{\overline{U}\Delta h''}$',...
            'C1:$\rho_a C_D\overline{\overline{U}T_o''\alpha \overline{\Delta h}}$',...
            'C2:$\rho_a C_D\overline{\overline{U}T_o''\alpha \Delta h''}$',...
            'D1:$\rho_a C_D\overline{T_o''\beta \overline{\Delta h}}$',...
            'D2:$\rho_a C_D\overline{T_o''\beta  \Delta h''}$',...
            'E1:$\rho_a C_D\overline{ (T_o'')^2\alpha\beta \overline{\Delta h}}$',...
            'E2:$\rho_a C_D\overline{ (T_o'')^2\alpha\beta \Delta h''}$'};
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
        
end

end










