plt_error_box = false;

t_range = 1:484;

year_vec = [2003];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

add_ssh_flag = false;

L = 250000 ; % m

alpha_pos_flag = false;

month_str = 'DJFM';

filter_type = 'fft'; % filter type 'lanczos' or 'boxcar'


%% begin
if alpha_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

for i = 1:length(year_vec)
    year=year_vec(i);
    load(sprintf('model_n_ERA_data_L_%d_%s%s_%s_%d',L/1000,con_str,month_str,filter_type,year))
    
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year));
    
    if plt_error_box
        
        err_box_lat = [32 38];
        err_box_lon = [140 160];
        
        err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
        err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));
        
        model_sshf = nanmean(model_full_sshf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        model_slhf = nanmean(model_full_slhf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        
        model_sshf_no_eddy = nanmean(model_no_eddy_sshf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        model_slhf_no_eddy = nanmean(model_no_eddy_slhf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        
        ERA5_sshf = nanmean(sshf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        ERA5_slhf = nanmean(slhf_patch(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        
        ERA5_sshf_CTRL = nanmean(era_no_eddy_sshf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        ERA5_slhf_CTRL = nanmean(era_no_eddy_slhf(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
        
        lat_er = lat(patch_lat);
        lat_er = lat_er(err_box_bnds_lat);
        lon_er = lon(patch_lon);
        lon_er = lon_er(err_box_bnds_lon);
        
    else
        
        
        
        if strcmp(filter_type,'fft')
            
            cf = (1/(2*L));
            box_opt = [32 38; 143 167];
            
            file_for_box = load(sprintf('Qs_QL_optimization_data_L_%d_filt_%s_%d',L/1000,filter_type,2003),'box_limits');
            
            box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
            box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);
            
            %             % to index out of *_patch fields
            opt_patch_lat = box_lat(patch_lat);
            opt_patch_lon = box_lon(patch_lon);
            
            prime_lat =  lat>=file_for_box.box_limits(1,1) & lat<=file_for_box.box_limits(1,2);
            prime_lon = lon>=file_for_box.box_limits(2,1) & lon<=file_for_box.box_limits(2,2);
            %
            %             % to index out of *_prime fields
            opt_prime_lat = box_lat(prime_lat);
            opt_prime_lon = box_lon(prime_lon);
            
            lat_er = lat(box_lat);
            lon_er = lon(box_lon);
            
            model_sshf = nanmean(model_full_sshf(:,:,t_range),3)';
            model_slhf = nanmean(model_full_slhf(:,:,t_range),3)';
            
            model_sshf_no_eddy = nanmean(model_no_eddy_sshf(:,:,t_range),3)';
            model_slhf_no_eddy = nanmean(model_no_eddy_slhf(:,:,t_range),3)';
            
            ERA5_sshf = nanmean(sshf_patch(opt_patch_lon,opt_patch_lat,t_range),3)';
            ERA5_slhf = nanmean(slhf_patch(opt_patch_lon,opt_patch_lat,t_range),3)';
            
            ERA5_sshf_CTRL = nanmean(era_no_eddy_sshf(:,:,t_range),3)';
            ERA5_slhf_CTRL = nanmean(era_no_eddy_slhf(:,:,t_range),3)';
            
        else
            model_sshf = nanmean(model_full_sshf(:,:,t_range),3)';
            model_slhf = nanmean(model_full_slhf(:,:,t_range),3)';
            
            model_sshf_no_eddy = nanmean(model_no_eddy_sshf(:,:,t_range),3)';
            model_slhf_no_eddy = nanmean(model_no_eddy_slhf(:,:,t_range),3)';
            
            ERA5_sshf = nanmean(sshf_patch(:,:,t_range),3)';
            ERA5_slhf = nanmean(slhf_patch(:,:,t_range),3)';
            
            ERA5_sshf_CTRL = nanmean(era_no_eddy_sshf(:,:,t_range),3)';
            ERA5_slhf_CTRL = nanmean(era_no_eddy_slhf(:,:,t_range),3)';
            
            lat_er = lat(patch_lat);
            lon_er = lon(patch_lon);
        end
        
    end
    
    
    % model_sshf = nanmean(sshf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    % model_slhf = nanmean(slhf_model(err_box_bnds_lon,err_box_bnds_lat,t_range),3)';
    %%
    ax = subplot(3,4,1);
    [~,h] = contourf(lon_er,lat_er,ERA5_sshf);
    title(sprintf('ERA5 full SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,2);
    [~,h] = contourf(lon_er,lat_er,ERA5_slhf);
    title(sprintf('ERA5 full SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,3);
    [~,h] = contourf(lon_er,lat_er,model_sshf);
    title(sprintf('model full SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,4);
    [~,h] = contourf(lon_er,lat_er,model_slhf);
    title(sprintf('model full SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,5);
    [~,h] = contourf(lon_er,lat_er,ERA5_sshf_CTRL);
    title(sprintf('ERA5 no eddy SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,6);
    [~,h] = contourf(lon_er,lat_er,ERA5_slhf_CTRL);
    title(sprintf('ERA5 no eddy SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,7);
    [~,h] = contourf(lon_er,lat_er,model_sshf_no_eddy);
    title(sprintf('no eddy SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,8);
    [~,h] = contourf(lon_er,lat_er,model_slhf_no_eddy);
    title(sprintf('no eddy SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(3,4,9);
    [~,h] = contourf(lon_er,lat_er,ERA5_sshf - ERA5_sshf_CTRL);
    max_diff = max(ERA5_sshf(:) - ERA5_sshf_CTRL(:));
    min_diff = min(ERA5_sshf(:) - ERA5_sshf_CTRL(:));
    title(sprintf('ERA5 diff SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax,max_diff,min_diff)
    
    ax = subplot(3,4,10);
    [~,h] = contourf(lon_er,lat_er,ERA5_slhf - ERA5_slhf_CTRL);
    max_diff = max(ERA5_slhf(:) - ERA5_slhf_CTRL(:));
    min_diff = min(ERA5_slhf(:) - ERA5_slhf_CTRL(:));
    title(sprintf('ERA5 diff SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax,max_diff,min_diff)
    
    ax = subplot(3,4,11);
    [~,h] = contourf(lon_er,lat_er,model_sshf-model_sshf_no_eddy);
    max_diff = max(model_sshf(:) - model_sshf_no_eddy(:));
    min_diff = min(model_sshf(:) - model_sshf_no_eddy(:));
    title(sprintf('model diff SSHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax,max_diff,min_diff)
    
    ax = subplot(3,4,12);
    [~,h] = contourf(lon_er,lat_er,model_slhf - model_slhf_no_eddy);
    max_diff = max(model_slhf(:) - model_slhf_no_eddy(:));
    min_diff = min(model_slhf(:) - model_slhf_no_eddy(:));
    title(sprintf('model diff SLHF $$[$$ Wm$$^{-2}]$$'),'interpreter','latex')
    format_fig(h,ax,max_diff,min_diff)
    
    set(gcf,'color','w','position',[ 232  1  1188  801],'NumberTitle','off','Name',num2str(year))
    
    if add_ssh_flag
        
        for i = 1:12
            subplot(3,4,i)
            hold on
            plot_SSH_contour;
        end
        
    end
    
    update_figure_paper_size()
    
    if add_ssh_flag
        print(sprintf('imgs/cmp_model_ERA5_L_%d_%sssh_%s_%d',L/1000,con_str,filter_type,year),'-dpdf')
    else
        print(sprintf('imgs/cmp_model_ERA5_L_%d_%s%s_%d',L/1000,con_str,filter_type, year),'-dpdf')
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












