
plt_error_box = false;

year = 2007;

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

add_ssh_flag = false;

month_str = 'DJFM';

alpha_pos_flag = false;

L = 500000;

%% begin

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

load(sprintf('flux_terms_%d_%s%d',L/1000,con_str,year))

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
    'lat','lon','patch_lat','patch_lon','time');

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

for j = 1:2
    tt_vec = false(length(time),1);
    if ismember('D',month_str)
        tt_vec = tt_vec | month(time)==12;
    end
    if ismember('J',month_str)
        tt_vec = tt_vec | month(time)==1;
    end
    if ismember('F',month_str)
        tt_vec = tt_vec | month(time)==2;
    end
    if ismember('M',month_str)
        tt_vec = tt_vec | month(time)==3;
    end
    p = sum(tt_vec);
    tt_inds = find(tt_vec);
    
    figure(j)
    ax = subplot(2,4,1);
    [~,h] = contourf(lon_er,lat_er,nanmean(term1(:,:,tt_inds,j),3)');
    title(sprintf('term1 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,2);
    [~,h] = contourf(lon_er,lat_er,nanmean(term2(:,:,tt_inds,j),3)');
    title(sprintf('term2 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,3);
    [~,h] = contourf(lon_er,lat_er,nanmean(term3(:,:,tt_inds,j),3)');
    title(sprintf('term3 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,4);
    [~,h] = contourf(lon_er,lat_er,nanmean(term4(:,:,tt_inds,j),3)');
    title(sprintf('term4 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,5);
    [~,h] = contourf(lon_er,lat_er,nanmean(term5(:,:,tt_inds,j),3)');
    title(sprintf('term5 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,6);
    [~,h] = contourf(lon_er,lat_er,nanmean(term6(:,:,tt_inds,j),3)');
    title(sprintf('term6 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,7);
    [~,h] = contourf(lon_er,lat_er,nanmean(term7(:,:,tt_inds,j),3)');
    title(sprintf('term7 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,8);
    [~,h] = contourf(lon_er,lat_er,nanmean(term8(:,:,tt_inds,j),3)');
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
        print(sprintf('imgs/flux_terms_%s_%s_ssh_L_%d_%d',title_str,month_str,L/1000,year),'-dpdf')
    else
        print(sprintf('imgs/flux_terms_%s_%s_L_%d_%d',title_str,month_str,L/1000,year),'-dpdf')
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






























