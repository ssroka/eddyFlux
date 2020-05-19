



plt_error_box = false;

t_range = 1:484;

year = 2003;

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

add_ssh_flag = false;

%% begin

load(sprintf('flux_terms_%d',year))


lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

for j = 1
    figure(j)
ax = subplot(2,4,1);
[~,h] = contourf(lon_er,lat_er,nanmean(term1(:,:,:,j),3)');
title(sprintf('term1 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,2);
[~,h] = contourf(lon_er,lat_er,nanmean(term2(:,:,:,j),3)');
title(sprintf('term2 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,3);
[~,h] = contourf(lon_er,lat_er,nanmean(term3(:,:,:,j),3)');
title(sprintf('term3 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,4);
[~,h] = contourf(lon_er,lat_er,nanmean(term4(:,:,:,j),3)');
title(sprintf('term4 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,5);
[~,h] = contourf(lon_er,lat_er,nanmean(term5(:,:,:,j),3)');
title(sprintf('term5 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,6);
[~,h] = contourf(lon_er,lat_er,nanmean(term6(:,:,:,j),3)');
title(sprintf('term6 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,7);
[~,h] = contourf(lon_er,lat_er,nanmean(term7(:,:,:,j),3)');
title(sprintf('term7 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

ax = subplot(2,4,8);
[~,h] = contourf(lon_er,lat_er,nanmean(term8(:,:,:,j),3)');
title(sprintf('term8 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
format_fig(h,ax)

if j == 1
    title_str = 'sensible';
else
    title_str = 'latent';
end

set(gcf,'color','w','position',[ 232  1  1188  801],...
    'NumberTitle','off','Name',sprintf('%s %d',title_str,year))



if add_ssh_flag
    
    for i = 1:8
        subplot(2,4,i)
        hold on
        plot_SSH_contour;
    end
    
end


update_figure_paper_size()

if add_ssh_flag
    print(sprintf('imgs/flux_terms_%s_ssh_%d',title_str,year),'-dpdf')
else
    print(sprintf('imgs/flux_terms_%s_%d',title_str,year),'-dpdf')
end

end

function [] = format_fig(h,plt_num,max_val,min_val)

set(h,'edgecolor','none')
set(gca,'ydir','normal','fontsize',15)
colorbar
xlabel('deg')
ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end



%{






%}






























