clear;close all;clc

filter_method = 'my_boxcar';
year = 2007;
term_num = 2;

filter_flag        = 'box'; % 'box' or 'zonal'
error_metric_flag  = 'point'; % 'box_sum' 'box_mean' or 'point'
err_box_lat = [32 38];
err_box_lon = [140 160];

patch_str = 'Kur'; % 'GS'   'Kur'

switch patch_str
    case 'GS'
        % Gulf Stream
        lat_bnds = [25 45];
        lon_bnds = [275 305];
    case 'Kur'
        % Kurishio
        lat_bnds = [25 45];
        lon_bnds = [130 170];
        lat_sm_bnds = [30 40];
end

% terms
titles  = {...
    'c_{pa} C_D^s U (T_o - T_a)'
    'c_{pa} C_D^s U''(T_o''-T_a'')'
    'c_{pa} C_D^s \alpha_s U T_o'' (T_o'' - T_a'')'
    'c_{pa} C_D^s \alpha_s U'' T_o'' (T_o-T_a)'
    'c_{pa} C_D^s \alpha_s U'' T_o'' (T_o'' - T_a'')'
    'c_{pa} C_D^s \alpha_s U T_o'' (T_o-T_a)'
    'c_{pa} C_D^s U'' (T_o - T_a)'
    'c_{pa} C_D^s U (T_o''-T_a'')'
    };
%%


load(sprintf('%d_model_eddy_min_no_eddy_results.mat',year),'patch_lon','patch_lat',...
    'patch_str','lat','lon','lat_bnds','lon_bnds','U_mag','DT_patch','Lv','c_p_air','rho_a',...
    'qo_patch','qa_patch','SST_prime','sshf_patch','slhf_eddy','sshf_eddy',...
    'slhf_patch','time','SST_patch','P0_patch','t2m_patch','RH_patch')
Nx = 19;
Ny = 19;

err_box_bnds_lat = (lat(patch_lat)>err_box_lat(1))&(lat(patch_lat)<err_box_lat(2));
err_box_bnds_lon = (lon(patch_lon)>err_box_lon(1))&(lon(patch_lon)<err_box_lon(2));

switch filter_method
    case 'my_boxcar'
        NaN_inds = isnan(t2m_patch(err_box_bnds_lon,err_box_bnds_lat,1));
        %         NaN_inds = false(159,79);
        [M] = my_smoother_mat(sum(err_box_bnds_lon),sum(err_box_bnds_lat),Ny,Nx,NaN_inds);
        filter_fxn =@(x) smooth_mat(x,M);
        
end

load(sprintf('X_FFINAL_%d_const_only',year),'X');
aCd0 = X{year-2002};
CD_ref_s = aCd0(3);
CD_ref_L = aCd0(4);

as_est =@(SST_prime) aCd0(1);
as = aCd0(1);
aL_est =@(SST_prime) aCd0(2);

nt = size(SST_patch,3);

lat_er = lat(patch_lat);
lat_er = lat_er(err_box_bnds_lat);
lon_er = lon(patch_lon);
lon_er = lon_er(err_box_bnds_lon);

U = U_mag(err_box_bnds_lon,err_box_bnds_lat,:);
U_bar = zeros(size(U));
U_prime = zeros(size(U));

term = zeros(size(U));
term_f = zeros(size(U));

T_o = SST_patch(err_box_bnds_lon,err_box_bnds_lat,:);
T_o_bar = zeros(size(T_o));
T_o_prime = zeros(size(T_o));

T_a = t2m_patch(err_box_bnds_lon,err_box_bnds_lat,:);
T_a_bar = zeros(size(T_a));
T_a_prime = zeros(size(T_a));

sum_of_all_terms = zeros(size(U,1),size(U,2));

for term_num = 1:8
    sum_of_term = zeros(size(U,1),size(U,2));
    for i = 1:5:nt
        T_o_bar(:,:,i) = filter_fxn(T_o(:,:,i));
        T_o_prime(:,:,i) = T_o(:,:,i)-T_o_bar(:,:,i);
        
        T_a_bar(:,:,i) = filter_fxn(T_a(:,:,i));
        T_a_prime(:,:,i) = T_a(:,:,i)-T_a_bar(:,:,i);
        
        U_bar(:,:,i) = filter_fxn(U(:,:,i));
        U_prime(:,:,i) = U(:,:,i)-U_bar(:,:,i);
        
        switch term_num
            case 1
                % c_pa C_Ds U (T_o - T_a)
                term(:,:,i) = c_p_air.*CD_ref_s.*U_bar(:,:,i).*(T_o_bar(:,:,i)-T_a_bar(:,:,i));
            case 2
                % c_pa C_Ds U'(T_o'-T_a')
                term(:,:,i) = c_p_air.*CD_ref_s.*U_prime(:,:,i).*(T_o_prime(:,:,i)-T_a_prime(:,:,i));
            case 3
                % c_pa C_Ds U alpha_s T_o' (T_o' - T_a') - cross terms make celsius or
                % Kelvin matter
                term(:,:,i) = c_p_air.*CD_ref_s.*U_bar(:,:,i).*as.*T_o_prime(:,:,i).*(T_o_prime(:,:,i)-T_a_prime(:,:,i));
            case 4
                % c_pa C_Ds U' alpha_s T_o' (T_o-T_a)
                term(:,:,i) = c_p_air.*CD_ref_s.*U_prime(:,:,i).*as.*T_o_prime(:,:,i).*(T_o_bar(:,:,i)-T_a_bar(:,:,i));
            case 5
                % c_pa C_Ds U' alpha_s T_o' (T_o' - T_a')
                term(:,:,i) = c_p_air.*CD_ref_s.*U_prime(:,:,i).*as.*T_o_prime(:,:,i).*(T_o_prime(:,:,i)-T_a_prime(:,:,i));
            case 6
                % c_pa C_Ds U alpha_s T_o' (T_o-T_a)
                term(:,:,i) = c_p_air.*CD_ref_s.*U_bar(:,:,i).*as.*T_o_prime(:,:,i).*(T_o_bar(:,:,i)-T_a_bar(:,:,i));
            case 7
                % c_pa C_Ds U' (T_o - T_a)
                term(:,:,i) = c_p_air.*CD_ref_s.*U_prime(:,:,i).*(T_o_bar(:,:,i)-T_a_bar(:,:,i));
            case 8
                % c_pa C_Ds U (T_o'-T_a')
                term(:,:,i) = c_p_air.*CD_ref_s.*U_bar(:,:,i).*(T_o_prime(:,:,i)-T_a_prime(:,:,i));
        end
        sum_of_term = sum_of_term + term(:,:,i);
        sum_of_all_terms = sum_of_all_terms + term(:,:,i);
    end
    
    
    
    
    subplot(3,4,term_num)
    contourf(lon_er,lat_er,filter_fxn(nanmean(term(:,:,:),3))')
    % xlabel('Deg. Lon')
    % ylabel('Deg. Lat')
    title(sprintf('$$%s$$ [W/m$$^2$$]',titles{term_num}),'interpreter','latex')
    set(gca,'ydir','normal','fontsize',20)
    colorbar
end

subplot(3,4,[9 10])
contourf(lon_er,lat_er,nanmean(sshf_eddy(err_box_bnds_lon,err_box_bnds_lat,:),3)')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
title('total model Q_s','interpreter','latex')
set(gca,'ydir','normal','fontsize',25)
colorbar

subplot(3,4,[11 12])
contourf(lon_er,lat_er,nanmean(sum_of_all_terms,3)')
xlabel('Deg. Lon')
ylabel('Deg. Lat')
title('term sum Q_s','interpreter','latex')
set(gca,'ydir','normal','fontsize',25)
colorbar

set(gcf,'position',[1          11        1425         793],'color','w')

% figure(2)
% 
% subplot(1,2,1)
% contourf(lon_er,lat_er,nanmean(term(:,:,:),3)')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% title(sprintf('term %d [W/m^2]',term_num))
% set(gca,'ydir','normal','fontsize',25)
% 
% colorbar
% subplot(1,2,2)
% contourf(lon_er,lat_er,nanmean(term_f(:,:,:),3)')
% xlabel('Deg. Lon')
% ylabel('Deg. Lat')
% title(sprintf('term %d filtered  [W/m^2]',term_num))
% colorbar
% 
% set(gca,'ydir','normal','fontsize',25)
% set(gcf,'color','w','position',[ -46          43        1450         762])
% set(gcf,'numbertitle','off','name',num2str(year))
% 
% 
% 
% 












