

tot_er = sshf_model+slhf_model-(sshf_ERA5_box+slhf_ERA5_box);
mean_abs_error = abs(nanmean(tot_er,3)./nanmean(sshf_ERA5_box+slhf_ERA5_box,3));
% std_abs_error = nanstd(abs(tot_er(:)./(sshf_ERA5_box(:)+slhf_ERA5_box(:))))

subplot(2,3,1+(i-1)*3)
contourf(lon_box,lat_box,mean_sshf_ERA5_box+mean_slhf_ERA5_box)
colorbar
x_tick = get(gca,'xtick');
x_tick_plot = cell(1,length(x_tick));
for j = 1:length(x_tick)
    x_tick_plot(j) = {sprintf('$$%d^{\\circ}$$',x_tick(i))};
end
y_tick = get(gca,'ytick');
y_tick_plot = cell(1,length(y_tick));
for j = 1:length(y_tick)
    y_tick_plot(j) = {sprintf('$$%d^{\\circ}$$',y_tick(i))};
end
set(gca,'ydir','normal','fontsize',15,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')

title(sprintf('$$\\langle Q^{ERA5}\\rangle$$ $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
xlabel('longitude','Interpreter','latex')
ylabel('latitude','Interpreter','latex')
clim_sshf = get(gca,'clim');

subplot(2,3,2+(i-1)*3)
contourf(lon_box,lat_box,mean_model_sshf+mean_model_slhf)
colorbar
set(gca,'ydir','normal','fontsize',15,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')
title(sprintf('$$\\langle Q^{%s}\\rangle$$ $$[$$ Wm$$^{-2}]$$',model_str_title),'interpreter','latex')
xlabel('longitude','Interpreter','latex')
ylabel('latitude','Interpreter','latex')

subplot(2,3,3+(i-1)*3)
contourf(lon_box,lat_box,mean_abs_error')
colorbar
set(gca,'ydir','normal','fontsize',15,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')
title(sprintf('$$\\frac{\\langle Q^{ERA5}\\rangle - \\langle Q^{%s}\\rangle}{\\langle Q^{ERA5}\\rangle}$$',model_str_title),'interpreter','latex')
xlabel('longitude','Interpreter','latex')
ylabel('latitude','Interpreter','latex')

set(gcf,'color','w','position',[1         200        1438         592],'NumberTitle','off','Name',num2str(year))

for plot_letters_i = 1:3
    subplot(2,3,plot_letters_i+(i-1)*3)
    th = text(0.03,0.1,[plot_letters(plot_letters_i) ')'],'units','normalized');
    set(th,'units','normalized','fontsize',20,'backgroundcolor','w')
end

update_figure_paper_size()
if all(abs(abCD_factor-1)<1e-10)
    print(sprintf('%simgs/cmp_model_ERA5_%d_%s_box%d_%s_%d',data_base,L/1000,filter_type,box_num,model_str,year),'-dpdf')
else
    print(sprintf('%simgs/cmp_model_ERA5_%d_%s_box%d_%s_%d_abCDFAC_%s',data_base,L/1000,filter_type,box_num,model_str,year,strrep(num2str(abCD_factor),'.','_')),'-dpdf')
end



















