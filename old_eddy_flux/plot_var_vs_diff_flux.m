

subplot(2,1,1)
point_name = sprintf('%d, L = %d km',M(ii_year,1),L/1000);
plot(M(ii_year,4),M(ii_year,2),s2p(ii_year),'markerfacecolor',c2p(ii_L),'displayname',point_name,'markersize',10)
title(sprintf('Mean DJFM Eddy - No Eddy'))
hold on
legend('-dynamiclegend')
set(gca,'ydir','normal','fontsize',20)
subplot(2,1,2)
plot(M(ii_year,4),M(ii_year,3),s2p(ii_year),'markerfacecolor',c2p(ii_L),'displayname',point_name,'markersize',10)
title(sprintf('Median DJFM Eddy - No Eddy'))
hold on
legend('-dynamiclegend')
set(gca,'ydir','normal','fontsize',20)
set(gcf,'color','w')
set(gcf,'color','w','position',[ 114         251        1317         547])




