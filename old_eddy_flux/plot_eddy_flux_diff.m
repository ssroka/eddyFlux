        subplot(2,3,1)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(sshf_eddy(:,inds_sm,:) + slhf_eddy(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Eddy Flux [%s]\n (mean DJFM %d)\n$$\\alpha$$ = % 4.3f $$\\qquad C_D^*$$ = %2.2e',error_units,time(end).Year,as_est(0),CD_ref);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        subplot(2,3,2)
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(sshf_sm(:,inds_sm,:) + slhf_sm(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('No Eddy Flux [%s]\n (mean DJFM %d)\n$$\\alpha$$ = % 4.3f $$\\qquad C_D^*$$ = %2.2e',error_units,time(end).Year,as_est(0),CD_ref);
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        subplot(2,3,3)
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(sshf_diff(:,inds_sm,:) + slhf_diff(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('Eddy minus no Eddy Flux');
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        subplot(2,3,4)
        inds_sm = find(patch_lat_sm(patch_lat));
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(sshf_diff(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('time avg SHF diff');
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        subplot(2,3,5)
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(slhf_diff(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('time avg LHF diff');
        title(t_str,'interpreter','latex')
        subplot(2,3,6)
        contourf(lon(patch_lon),lat(patch_lat_sm),nanmean(slhf_patch(:,inds_sm,:),3)'+nanmean(sshf_patch(:,inds_sm,:),3)',30,'k')
        colorbar
        xlabel(' Deg lon ')
        ylabel(' Deg lat ')
        p_str = '%';
        t_str = sprintf('total ERA flux');
        title(t_str,'interpreter','latex')
        set(gca,'ydir','normal','fontsize',20)
        set(gcf,'color','w')
        set(gcf,'color','w','position',[ 114         251        1317         547])

        