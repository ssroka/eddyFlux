% data from
%{
https://resources.marine.copernicus.eu/?option=com_csw&task=results
https://resources.marine.copernicus.eu/?option=com_csw&view=order&record_id=c0635fc4-07d3-4309-9d55-cfd3e6aa788b

%}
% clear;clc;close all

<<<<<<< HEAD

% data_base = '/Users/ssroka/MIT/Research/eddyFlux/';
data_base = '/Volumes/SydneySroka_Remy/eddyFlux/';
year_vec = [2003:2018];
L        = 250000; % m
filter_type = 'fft'; % filter type 'lanczos' or 'boxcar' or 'fft'
data_src = [data_base 'ERA5_data/'];
debug_flag = false;
% plot_flag = false;
calc_alpha_beta_CD_flag = true;
plot_model_ERA_err_flag = true;
abCD_factor = 1;
alpha_pos_flag = false;
beta_pos_flag = false;
fft_first_flag = true;
plot_all_alpha_beta_CD = false;
plot_all_ABC_vs_year = false;
plot_all_QandC = false;
plot_ABC_comp = false;

=======
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
box_num = 3;

switch box_num
    case 1
        box_opt = [36 41.5; 143 152];
    case 2
        box_opt = [30 44.5; 148 169];
    case 3
        box_opt = [30 41.5; 142.5 169];
    case 4
        box_opt = [30 41.5; 142.5 153];
    case 5
        box_opt = [30 41.5; 153 169];
end
<<<<<<< HEAD

model_str = 'beta';
=======
    
model_str = 'alpha';
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(1)),...
    'lat','lon','patch_lat','patch_lon','sshf_patch','slhf_patch');

% filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year);

% file_for_box = load(sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%d',L/1000,filter_type,box_num,year_vec(1)),'box_opt');
<<<<<<< HEAD
% filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year_vec(1));
% file_for_box = load(filename,'sshf_patch','slhf_patch');
=======
filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year_vec(1));
file_for_box = load(filename,'bs_multiplier','bL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch');
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
setup_lat_lon_vec;

box_lat = lat>=box_opt(1,1) & lat<=box_opt(1,2);
box_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

% to index out of *_patch fields
opt_patch_lat = box_lat(patch_lat);
opt_patch_lon = box_lon(patch_lon);

prime_lat =  lat>=box_opt(1,1) & lat<=box_opt(1,2);
prime_lon = lon>=box_opt(2,1) & lon<=box_opt(2,2);

% to index out of *_prime fields
opt_prime_lat = box_lat(prime_lat);
opt_prime_lon = box_lon(prime_lon);

lat_plot = lat(box_lat);
lon_plot = lon(box_lon);

m = sum(opt_prime_lon);
n = sum(opt_prime_lat);

if beta_pos_flag
    con_str = 'cons_';
else
    con_str = '';
end

% salinity = 34*ones(m,n);% ppt for Lv calculation

% if strcmp(filter_type,'boxcar')
%     [M] = boxcar_filter_mat(m,n,Ny,Nx,NaN_inds);
% end

<<<<<<< HEAD
=======
year_vec = [2003:2018];
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28


% data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';
addpath([data_base 'filter'])
addpath([data_base 'get_CD_alpha'])
addpath(['/Users/ssroka/Documents/MATLAB/util/'])
addpath(['/Users/ssroka/Documents/MATLAB/mpm/sandbox/NayarFxns'])

%     box_opt = [32 38; 143 167];
%     box_opt = [30 42; 144 168];
% box_opt = [30 44.5; 148 169];
% box_opt = [30 41.5; 142.5 169];
cntr = [0:0.2:0.8];
cntr_plot = [-0.2:0.2:1];

%% filter set up
load('env_const.mat')

d = zeros(length(cntr),length(year_vec));

for j = 1:length(year_vec)
    year = year_vec(j);
    
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year),...
        'SST_patch','lat','lon','patch_lat','patch_lon',...
        'slhf_patch','sshf_patch');
    
    cntr_length = zeros(length(cntr),length(year_vec));
    
    fprintf('year %d\n',year)
    
    % SSH_data = load('/Users/ssroka/MIT/Research/eddyFlux/CCAR_SSH_data/SSH_20000103_20090627_CCAR.mat','lon','lat','time','SSH');
    if year<2008 & year>2002
        SSH_data = load([data_base 'ssh_data_global/global-reanalysis-phy-001-030-daily_20021201_20070401.mat'],'lon','lat','time','SSH');
    elseif year<2013 & year>2007
        SSH_data = load([data_base 'ssh_data_global/global-reanalysis-phy-001-030-daily_20071201_20120331.mat'],'lon','lat','time','SSH');
    elseif year>2012 & year<2018
        SSH_data = load([data_base 'ssh_data_global/global-reanalysis-phy-001-030-daily_20131201_20180331.mat'],'lon','lat','time','SSH');
        
    end
    
    %     patch_lat_bnds = [min(lat(patch_lat)) max(lat(patch_lat))];
    %     patch_lon_bnds = [min(lon(patch_lon)) max(lon(patch_lon))];
    
    patch_lat_bnds = box_opt(1,1:2);
    patch_lon_bnds = box_opt(2,1:2);
    
    patch_lat_SSH = (SSH_data.lat>=patch_lat_bnds(1)) & (SSH_data.lat<=patch_lat_bnds(2));
    patch_lon_SSH = (SSH_data.lon>=patch_lon_bnds(1)) & (SSH_data.lon<=patch_lon_bnds(2));
    
    % switch year
    %     case 2003
    %         ssh_lvl = 0.5;
    %     case 2007
    %         ssh_lvl = 0.5;
    % end
    SSH_mat = zeros(sum(patch_lat_SSH),sum(patch_lon_SSH),length(SSH_data.time));
    count = 1;
    clear SSH_4_mean
    for i = 1:length(SSH_data.time)
        Dec_SSH = (SSH_data.time.Year(i)==year-1) && (SSH_data.time.Month(i) == 12);
        JFM_SSH = (SSH_data.time.Year(i)==year) && (SSH_data.time.Month(i) <= 3);
        if (Dec_SSH || JFM_SSH)
            SSH_mat(:,:,count) = SSH_data.SSH(patch_lon_SSH,patch_lat_SSH,i)';
            count = count + 1;
        end
    end
<<<<<<< HEAD
    SSH_mat = SSH_mat(:,:,1:count-1);
    %         meanSSH = mean(SSH,3);
    
    for ii_SSH = 1:size(SSH_mat,3)
        
        figure(year)
        [Cf,ch] = contourf(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
            SSH_mat(:,:,ii_SSH)); %this fills in the lower and upper regions
        hold on
        [Cf,ch] = contourf(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
            SSH_mat(:,:,ii_SSH),cntr_plot);
        clim = get(gca,'clim');
        hold on
        [C,h] = contour(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
            SSH_mat(:,:,ii_SSH),cntr,'w','linewidth',2);
        clabel(C,h,'color','w')
        set(h,'levellist',cntr)
        set(gca,'clim',clim)
        C(:,C(1,:)<20) = [];
        
        if ii_SSH==32 && j==1
            % set(gca,'edgecolor','none')
            set(gca,'ydir','normal','fontsize',15)
            colorbar
            set(gca,'fontsize',25)
            set(gcf,'color','w','position',[1 215 1439 587],'NumberTitle','off','Name',num2str(year))
            monthname = month(SSH_data.time(ii_SSH),'name');
            yearname = num2str(year);
            dayname = num2str(day(SSH_data.time(ii_SSH)));
            title(sprintf('Sea Surface Height Anomaly %s %s, %s [m]',monthname{1},dayname,yearname),'interpreter','latex')
            format_fig
            update_figure_paper_size()
            print(sprintf('imgs/ssh_cntr_%d',year),'-dpdf')
        end
        
        Vq = interp2(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
            SSH_mat(:,:,ii_SSH),C(1,:),C(2,:));
        for i = 1:length(cntr)
            %         bM = regionprops(meanSSH>=cntr(i),'Perimeter');
            %         cntr_length(i,j) = bM.Perimeter;
            inds = find(abs(Vq-cntr(i))<1e-10);
            for k = 1:length(inds)-1
                delta = norm(C(:,inds(k+1))-C(:,inds(k)));
                if delta < 0.25
                    d(i,j) = d(i,j) + delta;
                end
            end
        end
    end
    %     load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year),'sshf_patch','slhf_patch')
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(j)),...
        'sshf_patch','slhf_patch');
    
    SSHF_mean = nanmean(nanmean(nanmean(sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),3)));
    SLHF_mean = nanmean(nanmean(nanmean(slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),3)));
    
=======
    SSH = SSH(:,:,1:count-1);
    
    meanSSH = mean(SSH,3);
    
    figure(year)
    [~,ch] = contourf(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        meanSSH);
    clim = get(gca,'clim');
    hold on
    [C,h] = contour(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        meanSSH,cntr,'w','linewidth',2);
    clabel(C,h,'color','w')
    set(h,'levellist',cntr)
    set(gca,'clim',clim)
    colorbar
    set(gca,'fontsize',25)
    set(gcf,'color','w','position',[1 215 1439 587],'NumberTitle','off','Name',num2str(year))
    title(['SSH variance (colorbar), white contours = time mean',num2str(year)])
    
    C(:,C(1,:)<20) = [];
    Vq = interp2(SSH_data.lon(patch_lon_SSH),SSH_data.lat(patch_lat_SSH),...
        meanSSH,C(1,:),C(2,:));
    for i = 1:length(cntr)
        %         bM = regionprops(meanSSH>=cntr(i),'Perimeter');
        %         cntr_length(i,j) = bM.Perimeter;
        inds = find(abs(Vq-cntr(i))<1e-10);
        for k = 1:length(inds)-1
            delta = norm(C(:,k+1)-C(:,k));
            if delta < 0.25
                d(i,j) = d(i,j) + delta;
            end
        end
    end
%     load(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%s_%d',L/1000,con_str,fft_str,filter_type,box_num,model_str,year),'sshf_patch','slhf_patch')
    load(sprintf('%sERA5_patch_data_%d.mat',data_src,year_vec(j)),...
    'sshf_patch','slhf_patch');

    SSHF_mean = nanmean(nanmean(nanmean(sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),3)));
    SLHF_mean = nanmean(nanmean(nanmean(slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,:),3)));
    
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
    HF(j) =   SSHF_mean+  SLHF_mean;
    %{
        inds = C(1,:)<20;
        C(:,inds) = [];
        plot(C(1,:),C(2,:),'mo')
        
        
    %}
<<<<<<< HEAD
    
=======

>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
    
    %     update_figure_paper_size()
    %     print(sprintf('imgs/ssh_%d',year),'-dpdf')
end
<<<<<<< HEAD

d = d/size(SSH_mat,3);

save(sprintf('SSH_length_box%d_%d_%d',box_num,year_vec(1),year_vec(end)))

%%
%{
year_vec = [2003:2018];
box_num = 3;
load(sprintf('SSH_length_box%d_%d_%d',box_num,year_vec(1),year_vec(end)))
%}

function [] = format_fig()

x_tick = get(gca,'xtick');
x_tick_plot = cell(1,length(x_tick));
for i = 1:length(x_tick)
    x_tick_plot(i) = {sprintf('$$%d^{\\circ}$$',x_tick(i))};
end
y_tick = get(gca,'ytick');
y_tick_plot = cell(1,length(y_tick));
for i = 1:length(y_tick)
    y_tick_plot(i) = {sprintf('$$%d^{\\circ}$$',y_tick(i))};
end
set(gca,'ydir','normal','fontsize',25,...
    'xtick',x_tick,'XTickLabel',x_tick_plot,...
    'ytick',y_tick,'YTickLabel',y_tick_plot,...
    'TickLabelInterpreter','latex')
xlabel('longitude','interpreter','latex')
ylabel('latitude','interpreter','latex')
end
=======
d
figure
for i = 1:size(d,2)
    plot(cntr, d(:,i),'o-','linewidth',2,'displayname',num2str(year_vec(i)))
    hold on
end
legend('-dynamiclegend')
xlabel('contour height')
ylabel('distance [deg]')
set(gca,'fontsize',25,'xtick',cntr)
set(gcf,'color','w','position',[1 215 1439 587])


update_figure_paper_size()
print(sprintf('imgs/ssh_length_box%d_%d_%d',box_num,year_vec(1),year_vec(end)),'-dpdf')

%%
figure

plot(d(3,:),HF,'x','linewidth',2,'markersize',10)
hold on
for i = 1:length(year_vec)
   th = text(d(3,i)+0.25,HF(i),num2str(year_vec(i)));
   set(th,'fontsize',15)
end

% fit line
p = polyfit(d(3,:),HF,1);
h = plot(d(3,:),polyval(p,d(3,:)),'r','linewidth',2);
lh = legend(h,sprintf('m = %3.2f $W/m^2/deg$',p(1)));
set(lh,'interpreter','latex');
xlabel('SSH length of 0.4m contour [deg]','interpreter','latex')
ylabel('mean total flux $[W/m^2]$','interpreter','latex')
set(gca,'fontsize',25)
set(gcf,'color','w','position',[1 215 1439 587])

update_figure_paper_size()
print(sprintf('imgs/ssh_length_vs_HF_box%d_%s_%d_%d',box_num,model_str,year_vec(1),year_vec(end)),'-dpdf')
>>>>>>> 9c55f0c62a0fbd90d8c3a9a5f04c6e48ab356b28
