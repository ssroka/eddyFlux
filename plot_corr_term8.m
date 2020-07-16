
plt_error_box = false;

year_vec = [2003];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

addpath('/Users/ssroka/MIT/Research/eddyFlux/filter')
addpath('/Users/ssroka/MIT/Research/eddyFlux/filter/lanczosfilter/')
addpath('/Users/ssroka/MIT/Research/eddyFlux/get_CD_alpha')
addpath('/Users/ssroka/Documents/MATLAB/util/')
addpath('/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns')

add_ssh_flag = false;

month_str = 'DJFM';

alpha_pos_flag = false;

L = 500000;

filter_type = 'lanczos'; % filter type 'lanczos' or 'boxcar'

%% begin

if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'SST_patch','lat','lon','patch_lat','patch_lon','time');

m = size(SST_patch,1);
n = size(SST_patch,2);

d_lat = abs(lat(2)-lat(1));
d_lon = abs(lon(2)-lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

Nx = floor(L/dx)+mod(floor(L/dx),2)+1; % make Nx odd
Ny = floor(L/dy)+mod(floor(L/dx),2)+1; % make Ny odd

NaN_inds = isnan(SST_patch(:,:,1));

if strcmp(filter_type,'boxcar')
    [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
elseif strcmp(filter_type(1:7),'lanczos')
    cf = (1/(2*L));
end

for i = 1:length(year_vec)
    year = year_vec(i);
    
    load(sprintf('flux_terms_%d_%sfilt_%s_%d',L/1000,con_str,filter_type,year))
    
    load(sprintf('model_n_ERA_data_L_%d_%s%s_%s_%d',L/1000,con_str,month_str,filter_type,year))
    
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
        
        if j == 1
            title_str = 'sensible';
            model_flux = model_full_sshf;
        else
            title_str = 'latent';
            model_flux = model_full_slhf;
        end
        
        figure(j)
        ax = subplot(2,4,1);
        [~,h] = contourf(lon_er,lat_er,nanmean(model_flux(:,:,tt_inds),3)');
        title(sprintf('model flux $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        ax = subplot(2,4,2);
        [~,h] = contourf(lon_er,lat_er,nanmean(term1(:,:,tt_inds,j),3)');
        title(sprintf('term1 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        ax = subplot(2,4,3);
        [~,h] = contourf(lon_er,lat_er,nanmean(model_flux(:,:,tt_inds)-term1(:,:,tt_inds,j),3)');
        title(sprintf('model flux - term1 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        ax = subplot(2,4,4);
        [~,h] = contourf(lon_er,lat_er,nanmean(term8(:,:,tt_inds,j),3)');
        title(sprintf('term8 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        ax = subplot(2,4,5);
        model_flux_mean = nanmean(model_flux(:,:,tt_inds),3)';
        term1_mean = nanmean(term1(:,:,tt_inds,j),3)';
        term8_mean = nanmean(term8(:,:,tt_inds,j),3)';
        plot(model_flux_mean(:)-term1_mean(:),term8_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term8_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 8')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)

        ax = subplot(2,4,6);
        sum_terms = term1(:,:,tt_inds,j)+...
                    term2(:,:,tt_inds,j)+...
                    term3(:,:,tt_inds,j)+...
                    term4(:,:,tt_inds,j)+...
                    term5(:,:,tt_inds,j)+...
                    term6(:,:,tt_inds,j)+...
                    term7(:,:,tt_inds,j)+...
                    term8(:,:,tt_inds,j);
        [~,h] = contourf(lon_er,lat_er,nanmean(sum_terms,3)');
        title(sprintf('$$\\Sigma$$ terms $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        ax = subplot(2,4,7);
        [~,h] = contourf(lon_er,lat_er,nanmean(model_flux(:,:,tt_inds) - sum_terms,3)');
        title(sprintf('model - $$\\Sigma$$ terms $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        

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
            print(sprintf('imgs/corr_term8_%s_%s_ssh_L_%d_%s_%d',title_str,month_str,L/1000,filter_type,year),'-dpdf')
        else
            print(sprintf('imgs/corr_term8_%s_%s_L_%d_%s_%d',title_str,month_str,L/1000,filter_type,year),'-dpdf')
        end
        
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






























