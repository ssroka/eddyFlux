

for i = 1:length(year_vec)
    year = year_vec(i);
    
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile,'qo_patch','qa_patch','SST_patch','DT_patch')
    
    switch model_str
        case 'alpha'
            load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
               
            figure(1)
            ax = subplot(2,3,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_CTRL,3)');
            title('$\overline{U}$ $[m/s]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_prime,3)');
            title('$U''$ $[m/s]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(dh_CTRL,3)');
            title('$\overline{\Delta h}$ $[W/m^2]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,5);
            [~,h] = contourf(lon_box,lat_box,nanmean(dh_prime,3)');
            title('$\Delta h''$ $[W/m^2]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(SST,3)');
            title('$T_o$ $[K]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,6);
            [~,h] = contourf(lon_box,lat_box,nanmean(To_prime,3)');
            title('$T_o''$ $[K]$','interpreter','latex')
            format_fig(h,ax)
            
            set(gcf,'color','w','position',[61 221 1359 581])
            update_figure_paper_size()
            print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_comp_L_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
            
        case 'beta'
            
            
            load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
            
             
            figure(1)
            ax = subplot(2,3,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_bar,3)');
            title('$\overline{U}$ $[m/s]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_prime,3)');
            title('$U''$ $[m/s]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(dh_bar,3)');
            title('$\overline{\Delta h}$ $[W/m^2]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,5);
            [~,h] = contourf(lon_box,lat_box,nanmean(dh_prime,3)');
            title('$\Delta h''$ $[W/m^2]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(SST,3)');
            title('$T_o$ $[K]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,6);
            [~,h] = contourf(lon_box,lat_box,nanmean(To_prime,3)');
            title('$T_o''$ $[K]$','interpreter','latex')
            format_fig(h,ax)
            
            set(gcf,'color','w','position',[61 221 1359 581])
            update_figure_paper_size()
            print(sprintf('%simgs/ABC_comp_L_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
            %--------------------
            
            
            ax = figure(2);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_bar,3)');
            title(sprintf('$\\langle\\overline{U}\\rangle$ [m/s], DJFM %d',year),'interpreter','latex')
            format_fig(h,ax)
%             th = text(0.02,0.07,'a)','units','normalized');
%             set(th,'units','normalized','fontsize',20,'backgroundcolor','w')
            set(gcf,'color','w','position',[1     1   720   365])
            update_figure_paper_size()
            print(sprintf('%simgs/ABC_comp_U_L_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
            ax = figure(3);
            cntr_lvls = linspace(min(nanmean(To_prime,3),[],'all'),max(nanmean(To_prime,3),[],'all'),30);
            [~,h] = contourf(lon_box,lat_box,nanmean(To_prime,3)',cntr_lvls);
            title(sprintf('$\\langle T_o'' \\rangle$ [K], DJFM %d',year),'interpreter','latex')
            format_fig(h,ax)
%             th = text(0.02,0.07,'b)','units','normalized');
%             set(th,'units','normalized','fontsize',20,'backgroundcolor','w')
            set(gcf,'color','w','position',[1     1   720   365])
            update_figure_paper_size()
            print(sprintf('%simgs/ABC_comp_To_L_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
            
        case 'alphabeta'
            
            load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
             
            figure(1)
            ax = subplot(2,3,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_bar,3)');
            title('$\overline{U}$ $[m/s]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(U_prime,3)');
            title('$U''$ $[m/s]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(dh_bar,3)');
            title('$\overline{\Delta h}$ $[W/m^2]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,5);
            [~,h] = contourf(lon_box,lat_box,nanmean(dh_prime,3)');
            title('$\Delta h''$ $[W/m^2]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(SST,3)');
            title('$T_o$ $[K]$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,3,6);
            [~,h] = contourf(lon_box,lat_box,nanmean(To_prime,3)');
            title('$T_o''$ $[K]$','interpreter','latex')
            format_fig(h,ax)
            
            set(gcf,'color','w','position',[61 221 1359 581])
            update_figure_paper_size()
            print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_comp_L_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
    end
    
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

% set(h,'edgecolor','none')
% set(gca,'ydir','normal','fontsize',15)
% colorbar
% xlabel('deg')
% ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(gca,rwb_map([max_val 0 min_val],100))
else
    colormap(parula)
end
drawnow


end