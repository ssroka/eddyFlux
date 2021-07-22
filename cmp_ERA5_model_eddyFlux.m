

tot_er = sshf_model+slhf_model-(sshf_ERA5_box+slhf_ERA5_box);
mean_abs_error = abs(nanmean(tot_er,3)./nanmean(sshf_ERA5_box+slhf_ERA5_box,3));
% std_abs_error = nanstd(abs(tot_er(:)./(sshf_ERA5_box(:)+slhf_ERA5_box(:))))

subplot(1,3,1)
contourf(lon_box,lat_box,mean_sshf_ERA5_box+mean_slhf_ERA5_box)
colorbar
<<<<<<< HEAD
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

title(sprintf('$$\\left(Q^{ERA5}\\right)_{\\overline{t}}$$ $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('longitude','Interpreter','latex')
ylabel('latitude','Interpreter','latex')
clim_sshf = get(gca,'clim');
=======
set(gca,'ydir','normal','fontsize',15,'clim',clim_sshf)
title(sprintf('model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')


subplot(3,2,4)
contourf(lon_box,lat_box,mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15,'clim',clim_slhf)
title(sprintf('model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,5)
contourf(lon_box,lat_box,mean_sshf_ERA5_box-mean_model_sshf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 - model SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

subplot(3,2,6)
contourf(lon_box,lat_box,mean_slhf_ERA5_box - mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('ERA5 -model SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('deg')
ylabel('deg')

set(gcf,'color','w','position',[407  1  1013 801],'NumberTitle','off','Name',num2str(year))

update_figure_paper_size()
if abs(abCD_factor-1)<1e-10
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_%d_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year),'-dpdf')
else
    print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_%d_%s_box%d_%s_%d_abCDFAC_%s',L/1000,filter_type,box_num,model_str,year,strrep(num2str(abCD_factor),'.','_')),'-dpdf')
end
figure(1000)
subplot(2,2,1)
contourf(lon_box,lat_box,(mean_sshf_ERA5_box-mean_model_sshf)./mean_sshf_ERA5_box)
colorbar
set(gca,'ydir','normal','fontsize',15)
title(sprintf('(ERA5(x,y) - model(x,y))/ERA5(x,y) SSHF '),'interpreter','latex')
xlabel('deg')
ylabel('deg')
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28

subplot(1,3,2)
contourf(lon_box,lat_box,mean_model_sshf+mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')
title(sprintf('$$\\left(Q^{%s}\\right)_{\\overline{t}}$$ $$[$$ Wm$$^{-2}]$$',model_str_title),'interpreter','latex')
xlabel('longitude','Interpreter','latex')
ylabel('latitude','Interpreter','latex')

subplot(1,3,3)
contourf(lon_box,lat_box,mean_abs_error')
colorbar
set(gca,'ydir','normal','fontsize',15,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')
title(sprintf('$$\\frac{\\left(Q^{ERA5}\\right)_{\\overline{t}} - \\left(Q^{%s}\\right)_{\\overline{t}}}{\\left(Q^{%s}\\right)_{\\overline{t}}}$$',model_str_title,model_str_title),'interpreter','latex')
xlabel('longitude','Interpreter','latex')
ylabel('latitude','Interpreter','latex')

set(gcf,'color','w','position',[1         509        1438         296],'NumberTitle','off','Name',num2str(year))

update_figure_paper_size()
<<<<<<< HEAD
if all(abs(abCD_factor-1)<1e-10)
    print(sprintf('%simgs/cmp_model_ERA5_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
else
    print(sprintf('%simgs/cmp_model_ERA5_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,strrep(num2str(abCD_factor),'.','_')),'-dpdf')
end
=======
print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/cmp_model_ERA5_rel_er_%d_%s_box%d_%s_%d_abCDFAC_%s',L/1000,filter_type,box_num,model_str,year,strrep(num2str(abCD_factor),'.','_')),'-dpdf')


>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28



















