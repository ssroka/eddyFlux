clear; close all; clc

%{

Ran get_eddy_contribution for each of 2003 and 2007 except for one set used
the Nayar function for Lv and for the other used a constant value of
 2260000 for Lv.


%}

year_vec = [2003 2007];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';
load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
             'lat','lon','patch_lat','patch_lon');

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

for y = 1:length(year_vec)
    woLv = load(sprintf('model_n_ERA_data_no_Lv_fxn_%d.mat',year_vec(y)));
    wLv = load(sprintf('model_n_ERA_data_Lv_fxn_%d.mat',year_vec(y)));
    figure(y)
    ax = subplot(3,2,1);
    [~,h] = contourf(lon_er,lat_er,nanmean(woLv.model_full_slhf,3)');
    format_fig(h,ax)
    title('full model slhf, Lv = 2.26E6')
    ax = subplot(3,2,2);
    [~,h] = contourf(lon_er,lat_er,nanmean(woLv.model_no_eddy_slhf,3)');
    format_fig(h,ax)
    title('no eddy model slhf, Lv = 2.26E6')
    ax = subplot(3,2,3);
    [~,h] = contourf(lon_er,lat_er,nanmean(wLv.model_full_slhf,3)');
    format_fig(h,ax)
    title('full model slhf, Lv = Lv(SST)')
    ax = subplot(3,2,4);
    [~,h] = contourf(lon_er,lat_er,nanmean(wLv.model_no_eddy_slhf,3)');
    format_fig(h,ax)
    title('no eddy model slhf, Lv = Lv(SST)')
    ax = subplot(3,2,5);
    [~,h] = contourf(lon_er,lat_er,nanmean(woLv.model_full_slhf-wLv.model_full_slhf,3)');
    format_fig(h,ax)
    title('full model slhf, 2.26E6 - Lv(SST)')
    ax = subplot(3,2,6);
    [~,h] = contourf(lon_er,lat_er,nanmean(woLv.model_no_eddy_slhf-wLv.model_no_eddy_slhf,3)');
    format_fig(h,ax)
    title('no eddy model slhf, 2.26E6 - Lv')
    
    set(gcf,'color','w','position',[ 232  1  1188  801],'NumberTitle','off','Name',num2str(year_vec(y)))

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



