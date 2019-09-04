clear;close all;clc
%% user input

% Velocity Field options:
%    Chelton coupling coefficient = 'uT'
%    use constant U = u_bar       = 'constant'
%    use baroclinic wave          = 'baroclinic_wave'
velocity_field_flag = 'baroclinic_wave';

rand_seed = 0;

ctrl_flag           = 0; % 1 = no eddies, 0 = eddies

plot_flag           = 0;
plot_Heat_Flux      = 1;
plot_LHF_param_flag = 1;
Chelton_fig_4_flag  = 0;

% ---------- Parameters --------------
% Domain Parameters - from Foussard et al. (2019), table 2
L_km     = 9216; % km sq domain edge length
y_sst_km = 4500; % km
l_sst_km = 1000; % km
l_an_km  = 1500; % km
SST_bar  = 285;  % K
D_SST    = 20;   % K (total SST difference from North to South)
D_theta  = 10;   % K
P0  = 101300;    % Pa (surface pressure)

% wind-related
rho_a = 1.2;     % kg m^-3
c_p_air = 1000;  % J / kg / K
Lv = 2.26e6;     % J/kg

CD_ref = 1e-3;   % reference drag coefficient
u_bar = 4;       % m/s

% ocean-eddy related
eddy_size = [100 500]; % km [lower upper] bound
n_eddies  = 200;       % number of eddies to put in domain

% air-sea coupling related
uT_coeff = 0.2; % 0.2 - 0.44


%% create SST field
% ---------- Create Field --------------
N = 200;               % number of gridpoints
x = linspace(0,L_km,N);
[X, Y] = meshgrid(x,x);

[SST_field,SST_CTRL] = create_SST_field(ctrl_flag, SST_bar,D_SST,...
    y_sst_km,l_sst_km,l_an_km,L_km,eddy_size,X,Y,n_eddies,rand_seed);

T_prime = SST_field - SST_CTRL;
% T_prime = SST_field - repmat(mean(SST_field,2),1,size(X,2));
if plot_flag
    figure(1)
    contourf(x,x,T_prime,30)
    colorbar
    xlabel('km')
    ylabel('km')
    title('SST'' [K]')
    set(gca,'fontsize',18)
    set(gcf,'color','w')
    try
        colormap(othercolor('Mthermometercolors'))
    catch
        addpath('~/Documents/MATLAB/util/othercolor/')
        colormap(othercolor('Mthermometercolors'))
    end
end


%% create wind field
switch velocity_field_flag
    case 'uT'
%         U = u_bar/(L_km/2)*(Y-L_km/2)+uT_coeff*T_prime;
        U = u_bar+uT_coeff*T_prime;
        V = zeros(size(X));
        P = P0*ones(size(X));
    case 'constant'
        U = u_bar*ones(size(X));
        V = zeros(size(X));
        P = P0*ones(size(X));
    case 'baroclinic_wave'
        delta_P = 3000; % Pa
        centerline_phi = 45; % degrees N
        lambda    = 4000; % km baroclinic wavelength
        [U,V,P] = create_baroclinic_wave(P0,delta_P,centerline_phi,L_km,N,lambda,rho_a);
end
if plot_flag
    figure(20)
    imagesc(x,x,P./1000)% convert to kPa
    set(gca,'ydir','normal')
    hold on
    colormap(othercolor('Mthermometercolors'))
    d_quiver = 5;
    inds_2_plt = 1:d_quiver:N;
    qh = quiver(X(inds_2_plt,inds_2_plt),Y(inds_2_plt,inds_2_plt),...
        U(inds_2_plt,inds_2_plt)-u_bar,V(inds_2_plt,inds_2_plt));
    set(qh,'color','k')
    set(gca,'ydir','normal')
    xlabel('km')
    ylabel('km')
    title('Pressure [kPa], Velocity Anomaly ($$\vec{u}-\bar{u}$$)','interpreter','latex')
    set(gca,'fontsize',18)
    set(gcf,'color','w')
    colorbar
end

%% create stress field
N_param = 3;
for i_param = 1:N_param
    DT = 0.5;       % K
    RH = 0.8;       % relative humidity
    % coupling parameter for drag to temperature field
    a = 0.1;
    
    switch i_param
        case 1
            % parameter to test (loops through)
            test_param = [0.001 0.01 0.05 0.1 0.5 1.0 5.0];
            test_param_str = '\alpha';
            test_param_var_str = 'a';
            cmap_mat = parula(length(test_param));
        case 2
            % parameter to test (loops through)
            test_param = .60:.05:.95;
            test_param_str = 'RH';
            test_param_var_str = 'RH';
            cmap_mat = parula(length(test_param));
        case 3
            % parameter to test (loops through)
            test_param = [0.1 0.25 0.5 0.75 1];
            test_param_str = '\Delta T';
            test_param_var_str = 'DT';
            cmap_mat = parula(length(test_param));
    end
    
    
    zonal_SHF_diff = zeros(N,length(test_param));
    zonal_LHF_diff = zeros(N,length(test_param));
    for ind_avgHeat = 1:20
        for i = 1:length(test_param)
            eval(sprintf('%s = test_param(i);',test_param_var_str))
            
            CD = (1+a*T_prime)*CD_ref;
            
            tau_XY(:,:,1) = CD.*rho_a.*U.*sqrt(U.^2+V.^2);
            tau_XY(:,:,2) = CD.*rho_a.*V.*sqrt(U.^2+V.^2);
            
            if plot_flag
                figure(21)
                imagesc(x,x,SST_field)
                hold on
                inds_2_plt = 1:d_quiver:N;
                qh = quiver(X(inds_2_plt,inds_2_plt),Y(inds_2_plt,inds_2_plt),...
                    tau_XY(inds_2_plt,inds_2_plt,1),tau_XY(inds_2_plt,inds_2_plt,2));
                set(qh,'color','k')
                colorbar
                set(gca,'fontsize',18)
                title('Temperature [K]')
                set(gca,'fontsize',18)
                xlabel('km')
                ylabel('km')
                set(gcf,'color','w')
            end
            
            %% calculate curl/div of stress
            % compute curl in units of stress/m
            [cz,cav] = curl(X*1000,Y*1000,tau_XY(:,:,1),tau_XY(:,:,2));
            
            % compute divergence in units of sress/m
            d = divergence(X*1000,Y*1000,tau_XY(:,:,1),tau_XY(:,:,2));
            
            if plot_flag
                figure(3)
                subplot(1,2,1)
                imagesc(x,x,cz); colorbar
                hold on
                d_quiver = 5;
                inds_2_plt = 1:d_quiver:N;
                qh = quiver(X(inds_2_plt,inds_2_plt),Y(inds_2_plt,inds_2_plt),...
                    U(inds_2_plt,inds_2_plt),V(inds_2_plt,inds_2_plt));
                title('\nabla x \tau')
                set(gca,'fontsize',18)
                xlabel('km')
                ylabel('km')
                
                subplot(1,2,2)
                imagesc(x,x,d); colorbar
                hold on
                d_quiver = 5;
                inds_2_plt = 1:d_quiver:N;
                qh = quiver(X(inds_2_plt,inds_2_plt),Y(inds_2_plt,inds_2_plt),...
                    U(inds_2_plt,inds_2_plt),V(inds_2_plt,inds_2_plt));
                title('\nabla \cdot \tau')
                set(gca,'fontsize',18)
                xlabel('km')
                ylabel('km')
                set(gcf,'color','w')
            end
            
            %% calculate the SST gradient
            
            % to match Chelton Fig 4, take gradient in units of m
            [gradSSTX,gradSSTY] = gradient(T_prime,(x(2)-x(1))*1000);
            
            u_n = U./(sqrt(U.^2+V.^2));
            v_n = V./(sqrt(U.^2+V.^2));
            
            % find along-wind component
            aw_gradT_x = (gradSSTX.*u_n+gradSSTY.*v_n).*u_n;
            aw_gradT_y = (gradSSTX.*u_n+gradSSTY.*v_n).*v_n;
            
            along_wind = sqrt(aw_gradT_x.^2+aw_gradT_y.^2);
            
            % find cross-wind component
            cw_gradT_x = gradSSTX-(gradSSTX.*u_n+gradSSTY.*v_n).*u_n;
            cw_gradT_y = gradSSTY-(gradSSTX.*u_n+gradSSTY.*v_n).*v_n;
            
            cross_wind = sqrt(cw_gradT_x.^2+cw_gradT_y.^2);
            
            if plot_flag
                figure(4)
                subplot(1,2,1)
                imagesc(x,x,along_wind)
                colorbar
                title('along wind gradient [K/km]')
                set(gca,'fontsize',18)
                xlabel('km')
                ylabel('km')
                
                set(gcf,'color','w')
                
                subplot(1,2,2)
                imagesc(x,x,cross_wind)
                colorbar
                title('cross-wind gradient [K/km]')
                set(gca,'fontsize',18)
                xlabel('km')
                ylabel('km')
                
                set(gcf,'color','w')
            end
            
            %% Fig 4 from Chelton
            
            SST_ind = true(size(X))&((Y<1000)&(Y>0)&(X<5000)&(X>1500));
            
            cw_2_plt = cross_wind(SST_ind);
            cz_2_plt = cz(SST_ind);
            
            % fit a line
            if Chelton_fig_4_flag
                i_var = cw_2_plt(:)*(1000)*(100);
                d_var = abs(cz_2_plt(:))*1e7;
            else
                i_var = cw_2_plt(:);
                d_var = abs(cz_2_plt(:));
            end
            
            designX = [ones(length(i_var),1),i_var];
            b = designX\d_var;
            SSresid = sum((d_var-designX*b).^2);
            SStotal = (length(i_var)-1)*var(d_var);
            rsq = 1 - SSresid/SStotal;
            
            m_curl(i) = b(2);
            
            if plot_flag
                figure(5)
                subplot(1,2,1)
                if Chelton_fig_4_flag
                    plot(cw_2_plt(:)*(1000)*(100),abs(cz_2_plt(:))*1e7,'o')
                    xlabel('cross-wind SST gradient (deg C / 100 km)')
                    ylabel('curl (N m^-3 x 10^7)')
                else
                    plot(cw_2_plt(:),abs(cz_2_plt(:)),'o')
                    xlabel('cross-wind SST gradient (deg C / m)')
                    ylabel('curl (N m^-3)')
                end
                hold on
                set(gca,'fontsize',18)
                plot([min(i_var) max(i_var)],b(1)+b(2)*[min(i_var) max(i_var)],'r',...
                    'linewidth',2)
                title(sprintf('R^2 = %2.2f',rsq))
            end
            
            aw_2_plt = along_wind(SST_ind);
            d_2_plt = d(SST_ind);
            
            % fit a line
            if Chelton_fig_4_flag
                i_var = aw_2_plt(:)*(1000)*(100);
                d_var = abs(d_2_plt(:))*1e7;
            else
                i_var = aw_2_plt(:);
                d_var = abs(d_2_plt(:));
            end
            
            designX = [ones(length(i_var),1),i_var];
            b = designX\d_var;
            SSresid = sum((d_var-designX*b).^2);
            SStotal = (length(i_var)-1)*var(d_var);
            rsq = 1 - SSresid/SStotal;
            
            m_div(i) = b(2);
            
            if plot_flag
                figure(5)
                subplot(1,2,2)
                if Chelton_fig_4_flag
                    plot(aw_2_plt(:)*(1000)*(100),abs(d_2_plt(:))*1e7,'o','markersize',10)
                    xlabel('down-wind SST gradient (deg C / 100 km)')
                    ylabel('divergence (N m^-3 x 10^7)')
                else
                    plot(aw_2_plt(:),abs(d_2_plt(:)),'o','markersize',10)
                    xlabel('down-wind SST gradient (deg C / m)')
                    ylabel('divergence (N m^-3)')
                end
                hold on
                plot([min(i_var) max(i_var)],b(1)+b(2)*[min(i_var) max(i_var)],'r',...
                    'linewidth',2)
                title(sprintf('R^2 = %2.3f',rsq))
                set(gca,'fontsize',18)
                set(gcf,'position',[36 227 1230 578],'color','w')
            end
            
            % calculate heat fluxes
            Q_s = zeros(size(X,1),size(X,2),2);
            Q_L = zeros(size(X,1),size(X,2),2);
            
            for j = 1:2
                if j == 1
                    T = SST_CTRL;
                    str = 'no eddies';
                else
                    T = SST_field;
                    str = 'with eddies';
                end
                % ---- qo ----
                try
                    qo = SAM_qsatWater(T, P0) ;
                catch
                    addpath('~/Documents/MATLAB/util/')
                    qo = SAM_qsatWater(T, P0) ;
                end
                % ---- qa ----
                e_sat = SAM_psatWater(T-DT);
                e = RH*e_sat;
                r = 0.622 * e ./ max(e, P0-e);
                qa = r./(1+r);
                
                Q_s(:,:,j) = rho_a.*c_p_air.*CD.*sqrt(U.^2+V.^2).*(get_air_sea_DT(T,qo,qa,RH));
                Q_L(:,:,j) = rho_a.*Lv.*CD.*sqrt(U.^2+V.^2).*(qo-qa);
                
            end
            
            zonal_SHF_diff(:,i,ind_avgHeat) = mean((Q_s(:,:,2)-Q_s(:,:,1)),2);
            zonal_LHF_diff(:,i,ind_avgHeat) = mean((Q_L(:,:,2)-Q_L(:,:,1)),2);
            
            if (plot_flag || plot_Heat_Flux) && (ind_avgHeat==1)
                figure(7)
                
                subplot(N_param,2,1+2*(i_param-1))
                plot(mean((Q_s(:,:,2)-Q_s(:,:,1)),2),x,'linewidth',2,'color',cmap_mat(i,:),...
                    'displayname',sprintf('%s = %2.3f',test_param_str,test_param(i)))
                xlabel('W m^{-2}')
                ylabel('km')
                str1 = 'Zonal Mean $(Q_s^{eddy} - Q_s^{ctrl})$';
                title(str1,'interpreter','latex')
                lh = legend('-dynamiclegend');
                set(lh,'location','eastoutside')
                set(gca,'fontsize',18)
                hold on
                
                subplot(N_param,2,2+2*(i_param-1))
                plot(mean((Q_L(:,:,2)-Q_L(:,:,1)),2),x,'linewidth',2,'color',cmap_mat(i,:),...
                    'displayname',sprintf('%s = %2.3f',test_param_str,test_param(i)))
                xlabel('W m^{-2}')
                ylabel('km')
                str1 = 'Zonal Mean $(Q_L^{eddy} - Q_L^{ctrl})$';
                title(str1,'interpreter','latex')
                lh = legend('-dynamiclegend');
                set(lh,'location','eastoutside')
                set(gca,'fontsize',18)
                hold on
                set(gcf,'position',[237  9  1103  792],'color','w')
            end
        end
    end
    mean_zonal_SHF_diff = mean(zonal_LHF_diff,3);
    integrate_SHF = zeros(size(mean_zonal_SHF_diff,2),1);
    for ic = 1:size(mean_zonal_SHF_diff,2)
        integrate_SHF(ic) = trapz(x,mean_zonal_SHF_diff(:,ic)');
    end
    
    mean_zonal_LHF_diff = mean(zonal_LHF_diff,3);
    integrate_LHF = zeros(size(mean_zonal_LHF_diff,2),1);
    for ic = 1:size(mean_zonal_LHF_diff,2)
        integrate_LHF(ic) = trapz(x,mean_zonal_LHF_diff(:,ic)');
    end
    
    if plot_flag || plot_LHF_param_flag
        figure(8)
        subplot(1,2,1)
        plot([1:length(test_param)]/length(test_param),integrate_SHF,'o-',...
            'displayname',sprintf('%s, [%2.2f %2.2f]',test_param_str,min(test_param),max(test_param)))
        hold on
        legend('-dynamiclegend')
        set(gca,'xtick',[])
        title('$\int \Delta Q_s dy$ [W/m]','interpreter','latex')
        set(gca,'fontsize',18)
        
        subplot(1,2,2)
        plot([1:length(test_param)]/length(test_param),integrate_LHF,'o-',...
            'displayname',sprintf('%s, [%2.2f %2.2f]',test_param_str,min(test_param),max(test_param)))
        hold on
        legend('-dynamiclegend')
        set(gca,'xtick',[])
        title('$\int \Delta Q_L dy$ [W/m]','interpreter','latex')
        set(gca,'fontsize',18)
    end
end

%%

if plot_flag
    f = fit(test_param,abs(m_curl)','exp1');
    figure(6)
    semilogx(test_param,abs(m_curl),'o-','linewidth',2)
    hold on
    semilogx(test_param,abs(m_div),'x')
    ylabel('coupling coeff')
    legend('curl and cross-wind','div. and along-wind')
    title(sprintf('exp fit: %2.2f exp ( %2.2f $$%s$$)',f.a,f.b,test_param_str),'interpreter','latex')
    set(gca,'fontsize',18)
end
%
function [air_sea_DT] = get_air_sea_DT(T,qo,qa,RH)

% poly_fit_q_star_of_T
p    = [
    0.000000007253132
    -0.000008027905314
    0.003346650250875
    -0.622425571112314
    43.554695793466863];
dq_star_dT = 4*T.^3*p(1) + 3*T.^2*p(2) + 2*T.^1*p(3) + p(4);

air_sea_DT = (((qo-qa)-qo.*(1-RH))./RH)./dq_star_dT;

end



%{
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/SST_prime','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/SST','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/P_u_constant','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/P_u_uT','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/P_u_baroclinicWave','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/slope_stress_vs_SSTgrad_constant','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/slope_stress_vs_SSTgrad_uT','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/slope_stress_vs_SSTgrad_bw','-dpng')

print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/zonal_flux','-dpng')
print('/Users/ssroka/MIT/Research/eddyFlux/imgs_channel/flux_dep_RH_a_DT','-dpng')




%}