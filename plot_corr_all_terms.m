
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

filter_type = 'boxcar'; % filter type 'lanczos' or 'boxcar'

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
        ax = subplot(3,4,1);
        [~,h] = contourf(lon_er,lat_er,nanmean(model_flux(:,:,tt_inds),3)');
        title(sprintf('model flux $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
                
        ax = subplot(3,4,2);
        [~,h] = contourf(lon_er,lat_er,nanmean(term1(:,:,tt_inds,j),3)');
        title(sprintf('term1 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        ax = subplot(3,4,3);
        [~,h] = contourf(lon_er,lat_er,nanmean(model_flux(:,:,tt_inds)-term1(:,:,tt_inds,j),3)');
        title(sprintf('model flux - term1 $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
        format_fig(h,ax)
        
        model_flux_mean = nanmean(model_flux(:,:,tt_inds),3)';
        term1_mean = nanmean(term1(:,:,tt_inds,j),3)';
        term2_mean = nanmean(term2(:,:,tt_inds,j),3)';
        term3_mean = nanmean(term3(:,:,tt_inds,j),3)';
        term4_mean = nanmean(term4(:,:,tt_inds,j),3)';
        term5_mean = nanmean(term5(:,:,tt_inds,j),3)';
        term6_mean = nanmean(term6(:,:,tt_inds,j),3)';
        term7_mean = nanmean(term7(:,:,tt_inds,j),3)';
        term8_mean = nanmean(term8(:,:,tt_inds,j),3)';
        
        
        ax = subplot(3,4,5);
        plot(model_flux_mean(:)-term1_mean(:),term1_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term1_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 1')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)
        
        ax = subplot(3,4,6);
        plot(model_flux_mean(:)-term1_mean(:),term2_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term2_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 2')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)
        
        ax = subplot(3,4,7);
        plot(model_flux_mean(:)-term1_mean(:),term3_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term3_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 3')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)
        
        ax = subplot(3,4,8);
        plot(model_flux_mean(:)-term1_mean(:),term4_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term4_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 4')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)

        ax = subplot(3,4,9);
        plot(model_flux_mean(:)-term1_mean(:),term5_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term5_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 5')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)
        
        ax = subplot(3,4,10);
        plot(model_flux_mean(:)-term1_mean(:),term6_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term6_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 6')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)
        
        ax = subplot(3,4,11);
        plot(model_flux_mean(:)-term1_mean(:),term7_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term7_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 7')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)
        
        ax = subplot(3,4,12);
        plot(model_flux_mean(:)-term1_mean(:),term8_mean(:),'o');
        [B,BINT,R,RINT,STATS] = regress(term8_mean(:),[ones(size(term8_mean(:))) model_flux_mean(:)-term1_mean(:)]);
        xlabel('full - term 1')
        ylabel('term 8')
        title(sprintf('$$R^2$$ statistic = %f',STATS(1)),'interpreter','latex')
        set(gca,'ydir','normal','fontsize',15)        

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
            print(sprintf('imgs/flux_term_correlation_%s_%s_ssh_L_%d_%s_%d',title_str,month_str,L/1000,filter_type,year),'-dpdf')
        else
            print(sprintf('imgs/flux_term_correlation_%s_%s_L_%d_%s_%d',title_str,month_str,L/1000,filter_type,year),'-dpdf')
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






























