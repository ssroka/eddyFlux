

plot_set_names = {'Ta_{LP}','To_{LP}','Ta_{HP}','To_{HP}'};

for i = 1:2
    subplot(1,2,i)
    scatter(plot_set_vec(:,1+(i-1)*2),plot_set_vec(:,2+(i-1)*2))
    xlabel(plot_set_names{1+(i-1)*2})
    ylabel(plot_set_names{2+(i-1)*2})
    title(['$$' plot_set_names{i} '$$   ' datestr(time(tt),'dd-mm-yyyy HH:MM:SS')],'interpreter','latex')
    colorbar
    set(gca,'ydir','normal','fontsize',20)
end
set(gcf,'color','w','position',[ 84          50        1156         7329])
