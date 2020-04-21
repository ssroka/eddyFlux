

plot_set_names = {'Ta_{LP}','To_{LP}','Ta_{HP}','To_{HP}'};

for i = 1:4
    subplot(2,2,i)
    contourf(lon(patch_lon),lat(patch_lat_sm),plot_set(:,inds_sm,i)',30,'k')
    title(['$$' plot_set_names{i} '$$   ' datestr(time(tt),'dd-mm-yyyy HH:MM:SS')],'interpreter','latex')
    colorbar
    xlabel(' Deg lon ')
    ylabel(' Deg lat ')
    set(gca,'ydir','normal','fontsize',20)
end
set(gcf,'color','w','position',[ 84          50        1156         7329])
