

for i = 1:length(year_vec)
    year = year_vec(i);
    
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile,'qo_patch','qa_patch','SST_patch','DT_patch')
    
    switch model_str
        case 'alpha'
            load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
            
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
                set(gca,'clim',[-1 1])
            end
            set(gcf,'color','w','position',[61 221 1359 581])
            update_figure_paper_size()
            print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_L_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
            
        case 'beta'
            
            
            load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
            
            
            figure(1)
            ax = subplot(2,2,1);
            [~,h] = contourf(lon_box,lat_box,nanmean(Ab,3)');
            title('A:$\overline{C_D\overline{U}\hspace{2mm}\overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,2,2);
            [~,h] = contourf(lon_box,lat_box,nanmean(Bb,3)');
            title('B:$\overline{C_D\overline{U}\Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,2,3);
            [~,h] = contourf(lon_box,lat_box,nanmean(Cb,3)');
            title('C:$\overline{C_D\beta T_o'' \Delta h''}$','interpreter','latex')
            format_fig(h,ax)
            
            ax = subplot(2,2,4);
            [~,h] = contourf(lon_box,lat_box,nanmean(Db,3)');
            title('D:$\overline{C_D\beta T_o''\overline{\Delta h}}$','interpreter','latex')
            format_fig(h,ax)
            
            set(gcf,'color','w','position',[61 221 1359 581])
            update_figure_paper_size()
            print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/ABC_L_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
            
        case 'alphabeta'
            
            load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year))
            
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