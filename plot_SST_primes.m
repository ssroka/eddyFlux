  figure(1)
        subplot(2,3,4)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),SST_smooth(:,inds_sm,tt)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        title('No Eddy SST Field [K]')
        set(gca,'ydir','normal','fontsize',20)
        subplot(2,3,5)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm), SST_smooth_CTRL(:,inds_sm)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        title('Zonally Uniform SST (''no eddies'') [K]')
        set(gca,'ydir','normal','fontsize',20)
        subplot(2,3,6)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),SST_prime_smooth(:,inds_sm,tt)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        title('SST'' for ''no eddy'' case [K]')
        set(gca,'ydir','normal','fontsize',20)
        subplot(2,3,1)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),SST_patch(:,inds_sm,tt)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        title('Full SST Field [K]')
        set(gca,'ydir','normal','fontsize',20)
        subplot(2,3,2)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm), SST_patch(:,inds_sm,tt)'-SST_prime(:,inds_sm,tt)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        title('boxcar smoothed SST (eddy) [K]')
        set(gca,'ydir','normal','fontsize',20)
        subplot(2,3,3)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),SST_prime(:,inds_sm,tt)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        t_str = sprintf('SST'' for eddy case \nfull minus boxcar smoothed [K]');
        title(t_str);
        set(gca,'ydir','normal','fontsize',20)        
        set(gcf,'color','w','position',[114         376        1220         422])