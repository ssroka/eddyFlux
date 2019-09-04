
a_lim = 5; 
Nr = 4;
Nc = 2;


% as_cutoff  = as;
% as_cutoff(as>a_lim) = a_lim;
% as_cutoff(as<-a_lim) = -a_lim;

as_sm = smooth2a(nanmedian(as,3),Nr,Nc);

% aL_cutoff  = aL;
% aL_cutoff(aL>a_lim) = a_lim;
% aL_cutoff(aL<-a_lim) = -a_lim;

aL_sm = smooth2a(nanmedian(aL,3),Nr,Nc);



subplot(2,3,1)
imagesc(lon(patch_lon),lat(patch_lat),nanmean(as,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
set(gca,'fontsize',18)
t_str = sprintf('$$\\alpha_s$$ (mean DJFM %d)\n from full field ERA5',time(end).Year);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,2)
imagesc(lon(patch_lon),lat(patch_lat),nanmedian(as,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
set(gca,'fontsize',18)
t_str = sprintf('$$\\alpha_s$$ \n cut off at $$\\alpha_s=\\pm$$%2.2f',a_lim);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,3)
imagesc(lon(patch_lon),lat(patch_lat),nanmean(as_sm,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
set(gca,'fontsize',18)
t_str = sprintf('$$\\alpha_s$$ box-car smoothed\n boxcar lat = %2.1f deg, boxcar lon = %2.1f deg',Nc*0.25,Nr*0.25);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,4)
imagesc(lon(patch_lon),lat(patch_lat),nanmean(aL,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
set(gca,'fontsize',18)
t_str = sprintf('$$\\alpha_L$$ (mean DJFM %d)\n from full field ERA5',time(end).Year);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,5)
imagesc(lon(patch_lon),lat(patch_lat),nanmedian(aL,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
set(gca,'fontsize',18)
t_str = sprintf('$$\\alpha_L$$ \n cut off at $$\\alpha_L=\\pm$$%2.2f',a_lim);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)

subplot(2,3,6)
imagesc(lon(patch_lon),lat(patch_lat),nanmean(aL_sm,3)')
colorbar
xlabel(' Deg lon ')
ylabel(' Deg lat ')
set(gca,'fontsize',18)
t_str = sprintf('$$\\alpha_L$$ box-car smoothed\n boxcar lat = %2.1f deg, boxcar lon = %2.1f deg',Nc*0.25,Nr*0.25);
title(t_str,'interpreter','latex')
set(gca,'ydir','normal','fontsize',20)
 
set(gcf,'position',[30          68        1293         730],'color','w')

%{

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_smooth','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/alpha_smooth_07','-dpng')



%}
