
letter_vec = 'abcd';
for i = 1:length(year_vec)
    year = year_vec(i);
    load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%d_%s_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,year_vec(i)),'abCD');
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile,'qo_patch','qa_patch','SST_patch','DT_patch')
    load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%d_%s_%d%s',L/1000,con_str,fft_str,filter_type,box_num,er_box_num,model_str,year,abCD_fac_str))
    %             load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
    
    switch model_str
        case 'alpha'
            
            figure(1)
            ax = subplot(2,4,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(A,3)');
            title('A:$\rho_a C_D\overline{U}\overline{\Delta h}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(B,3)');
            title('B:$\rho_a C_D\overline{U''\Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(C1,3)');
            title('C1:$\rho_a C_D\overline{U}\overline{T_o''\alpha \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(C2,3)');
            title('C2:$\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,5);
            [~,h] = contourf(lon_box,lat_box,nanmean(C3,3)');
            title('C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,6);
            [~,h] = contourf(lon_box,lat_box,nanmean(D,3)');
            title('D:$\rho_a C_D\overline{U''\overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,7);
            [~,h] = contourf(lon_box,lat_box,nanmean(E1,3)');
            title('E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,8);
            [~,h] = contourf(lon_box,lat_box,nanmean(E2,3)');
            title('E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            for i = 2:8
                subplot(2,4,i)
                %                 set(gca,'clim',[-1 1])
            end
            
            
        case 'beta'
            ax = figure(2);
            dyn_thermo_ratio = abCD(1).*nanmean(To_prime,3)'./nanmean(U_bar,3)';
            [~,h] = contourf(lon_box,lat_box,dyn_thermo_ratio,40);
            %             [~,h] = contourf(lon_box,lat_box,(nanmean(Cb,3)'+nanmean(Db,3)')./(nanmean(Ab,3)'+nanmean(Bb,3)'),40);
            title(sprintf('$\\frac{\\langle\\beta T_o''\\rangle}{\\langle\\overline{U}\\rangle}$, DJFM %d',year),'interpreter','latex')
            format_fig(h,ax)
            set(gca,'ydir','normal','fontsize',20)
            %             th = text(0.02,0.07,'a)','units','normalized');
            %             set(th,'units','normalized','fontsize',20,'backgroundcolor','w')
            set(gcf,'color','w','position',[1     1   720   365])
            update_figure_paper_size()
            print(sprintf('%simgs/ABC_ratio_L_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,abCD_fac_str),'-dpdf')
            
            
            figure(3)
            [~,h] = contourf(lon_box,lat_box,nanmean(Db+Cb,3)');
            title(sprintf('$\\langle \\overline{Q^\\beta_3} + \\overline{Q^\\beta_4}\\rangle $ $$[$$ Wm$$^{-2}]$$, DJFM %d',year),'interpreter','latex')
            format_fig(h,ax)
            set(gcf,'color','w','position',[1     1   720   365])
            update_figure_paper_size()
            print(sprintf('%simgs/ABC_dynamic_L_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,abCD_fac_str),'-dpdf')
            
            
            
            figure(1)
            ax = subplot(2,2,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(Ab,3)');
            title('$\langle \overline{Q^\beta_1}\rangle $: $\overline{\rho_a C_D\overline{U}\hspace{1mm}\overline{\Delta h}}$ $$[$$ Wm$$^{-2}]$$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,2,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(Bb,3)');
            title('$\langle \overline{Q^\beta_2}\rangle $: $\overline{\rho_a C_D\overline{U}\Delta h''}$ $$[$$ Wm$$^{-2}]$$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,2,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(Cb,3)');
            title('$\langle \overline{Q^\beta_4}\rangle $: $\overline{\rho_a C_D\beta T_o'' \Delta h''}$ $$[$$ Wm$$^{-2}]$$','interpreter','latex')
            format_fig(h,ax)
            set(gca,'clim',[-1.5 2.5])
            
            ax = subplot(2,2,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(Db,3)');
            title('$\langle \overline{Q^\beta_3}\rangle $: $\overline{\rho_a C_D\beta T_o''\overline{\Delta h}}$ $$[$$ Wm$$^{-2}]$$','interpreter','latex')
            format_fig(h,ax)
            set(gca,'clim',[-1.5 2.5])
            
            
            for plot_num = 1:4
                ax = subplot(2,2,plot_num);
                th = text(0.07,0.15,sprintf('%s)',letter_vec(plot_num)),'units','normalized');
                set(th,'fontsize',20,'backgroundcolor','w')
            end
            
            %             sum_pos = get(gca,'position');
            %             vert_pos = 0.5-0.5*sum_pos(4);
            %             set(gca,'clim',[-2 4],'position',[0.7226 vert_pos 0.1810 0.3381])
            
            
            
        case 'alphabeta'
            ax = figure(2);
            dyn_thermo_ratio =  (nanmean(U_bar.*To_prime,3)'*abCD(1)+...
                abCD(2).*nanmean(To_prime,3)'+...
                abCD(1).*abCD(2).*(nanmean(To_prime.^2,3)'))./nanmean(U_bar,3)';
            [~,h] = contourf(lon_box,lat_box,dyn_thermo_ratio,40);
            title(sprintf('$\\left(\\langle\\overline{U}\\alpha T_o''\\rangle+\\langle\\beta T_o''\\rangle + \\langle\\alpha\\beta \\left(T_o''\\right)^2\\rangle\\right)/\\langle\\overline{U}\\rangle$, DJFM %d',year),'interpreter','latex')
            format_fig(h,ax)
            
            set(gca,'ydir','normal','fontsize',20)
            set(gcf,'color','w','position',[1     1   720   365])
            update_figure_paper_size()
            print(sprintf('%simgs/ABC_ratio_L_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,abCD_fac_str),'-dpdf')
            
            
            figure(1)
            ax = subplot(2,4,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(Aab,3)');
            title('A:$\rho_a C_D\overline{\overline{U}\overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(Bab,3)');
            title('B:$\rho_a C_D\overline{\overline{U}\Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(C1ab,3)');
            title('C1:$\rho_a C_D\overline{\overline{U}T_o''\alpha \overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(C2ab,3)');
            title('C2:$\rho_a C_D\overline{\overline{U}T_o''\alpha \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,5);
            [~,h] = contourf(lon_box,lat_box,nanmean(D1ab,3)');
            title('D1:$\rho_a C_D\overline{T_o''\beta \overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,6);
            [~,h] = contourf(lon_box,lat_box,nanmean(D2ab,3)');
            title('D2:$\rho_a C_D\overline{T_o''\beta  \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,7);
            [~,h] = contourf(lon_box,lat_box,nanmean(E1ab,3)');
            title('E1:$\rho_a C_D\overline{ (T_o'')^2\alpha\beta \overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,4,8);
            [~,h] = contourf(lon_box,lat_box,nanmean(E2ab,3)');
            title('E2:$\rho_a C_D\overline{ (T_o'')^2\alpha\beta \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            for i = 2:8
                subplot(2,4,i)
                %                 set(gca,'clim',[-1 1])
            end
            
    end
    set(gcf,'color','w','position',[61 221 1359 581])
    update_figure_paper_size()
    print(sprintf('%simgs/ABC_L_%d_%s_box%d_%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,er_box_num,model_str,year,abCD_fac_str),'-dpdf')
end


function [] = format_fig(h,plt_num,max_val,min_val)

set(h,'edgecolor','none')
set(gca,'ydir','normal','fontsize',15)
colorbar
x_tick = get(gca,'xtick');
x_tick_plot = cell(1,length(x_tick));
for i = 1:length(x_tick)
    x_tick_plot(i) = {sprintf('$$%d^{\\circ}$$',x_tick(i))};
end
y_tick = get(gca,'ytick');
y_tick_plot = cell(1,length(y_tick));
for i = 1:length(y_tick)
    y_tick_plot(i) = {sprintf('$$%d^{\\circ}$$',y_tick(i))};
end
set(gca,'ydir','normal','fontsize',15,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')
xlabel('longitude','interpreter','latex')
ylabel('latitude','interpreter','latex')

if nargin>2 % red white and blue colormap
    colormap(gca,rwb_map([max_val 0 min_val],100))
else
    colormap(parula)
end
drawnow


end