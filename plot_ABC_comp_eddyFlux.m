

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
            print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_comp_L_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
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
xlabel('deg')
ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(gca,rwb_map([max_val 0 min_val],100))
else
    colormap(parula)
end
drawnow


end