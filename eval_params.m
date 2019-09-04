clear;close all;clc
%% user input

use_uT_coeff_flag = 0;
ctrl_flag = 0; % 1 = no eddies, 0 = eddies

%% create SST field

% ====================================
% ============= Eddies ===============
% ====================================

% ---------- Parameters --------------
% Foussard et al. (2019)
% table 2
L_km = 9216;     % km sq domain edge length
y_sst_km = 4500; % km
l_sst_km = 1000; % km
l_an_km  = 1500; % km
SST_bar  = 285;  % K
D_SST    = 20;   % K
D_theta  = 10;   % K
P0  = 101300;     % hPa

eddy_size = [100 500]; % km [lower upper] bound
n_eddies  = 200;       % number of eddies to put in domain

% ---------- Create Field --------------
N = 200;               % number of gridpoints
x = linspace(0,L_km,N);
[X, Y] = meshgrid(x,x);

[SST_field,SST_CTRL] = create_SST_field(ctrl_flag, SST_bar,D_SST,...
    y_sst_km,l_sst_km,l_an_km,L_km,eddy_size,X,Y,n_eddies);

T_prime = SST_field - SST_CTRL;
% T_prime = SST_field - repmat(mean(SST_field,2),1,size(X,2));

figure(1)
contourf(x,x,SST_field,30)
colorbar
xlabel('km')
ylabel('km')
title('SST [K]')
set(gca,'fontsize',18)
try
    colormap(othercolor('RdBu9'))
catch
    addpath('~/Documents/MATLAB/util/othercolor/')
    colormap(othercolor('RdBu9'))
end
%% create wind field
rho_a = 1.2;    % kg m^-3
u_bar = 4;    % m/s

if use_uT_coeff_flag
    uT_coeff = 0.2; % 0.2 - 0.44
    U = u_bar/(L_km/2)*(Y-L_km/2)+uT_coeff*T_prime;
    V = zeros(size(X));
else
    U = u_bar*ones(size(X));
    V = zeros(size(X));
end

%% create stress field
as = [0.001 0.01 0.05 0.1 0.25 0.5 1];
for i = 1:length(as)
    a = as(i);
    CD_ref = 1e-3;
    
    CD = (1+a*T_prime)*CD_ref;
    
    tau_XY(:,:,1) = CD.*rho_a.*U.*sqrt(U.^2+V.^2);
    tau_XY(:,:,2) = CD.*rho_a.*V.*sqrt(U.^2+V.^2);
    
    figure(2)
    subplot(1,2,1)
    imagesc(tau_XY(:,:,1))
    title('\tau_x [N m^-2]')
    colorbar
    set(gca,'fontsize',18)

    subplot(1,2,2)
    imagesc(tau_XY(:,:,1))
    colorbar
    title('\tau_y [N m^-2]')
    set(gca,'fontsize',18)

    
    %% calculate curl/div of stress
    figure(3)
    % compute curl
    [cz,cav] = curl(X*1000,Y*1000,tau_XY(:,:,1),tau_XY(:,:,2));
    subplot(1,2,1)
    imagesc(x,x,cz); colorbar
    hold on
    q = quiver(X,Y,U,V);
    title('\nabla x \tau')
    set(gca,'fontsize',18)
    
    
    d = divergence(X*1000,Y*1000,tau_XY(:,:,1),tau_XY(:,:,2));
    subplot(1,2,2)
    imagesc(x,x,d); colorbar
    hold on
    q = quiver(X,Y,U,V);
    title('\nabla \cdot \tau')
    set(gca,'fontsize',18)
    
    %% calculate the SST gradient
    
    [gradSSTX,gradSSTY] = gradient(SST_field,(x(2)-x(1))/100);
    
    figure(4)
    subplot(1,2,1)
    imagesc(x,x,gradSSTX)
    colorbar
    title('along wind gradient')
    set(gca,'fontsize',18)

    subplot(1,2,2)
    imagesc(x,x,gradSSTY)
    colorbar
    title('cross-wind gradient')
    set(gca,'fontsize',18)

    %% Fig 4 from Chelton
    
    figure(5)
    subplot(1,2,1)
    plot(gradSSTY(:),cz(:)*1e7,'o')
    xlabel('cross-wind SST gradient (deg C / 100 km)')
    ylabel('curl (N m^-3 x 10^7)')
    % calc slope
    p = polyfit(gradSSTY(:),cz(:)*1e7,1);
    m_curl(i) = p(1);% (cz(end)-cz(1))*1e7/(gradSSTY(end)-gradSSTY(1));
    title(['slope = ' num2str(m_curl(i))])
    set(gca,'fontsize',18)

    subplot(1,2,2)
    plot(gradSSTX(:),d(:)*1e7,'o')
    xlabel('along-wind SST gradient (deg C / 100 km)')
    ylabel('divergence (N m^-3 x 10^7)')
    m_div(i) = (d(end)-d(1))*1e7/(gradSSTX(end)-gradSSTX(1));
    title(['slope = ' num2str(m_div(i))])
    set(gca,'fontsize',18)

end
%%
figure(6)
semilogx(as,abs(m_curl),'o-','linewidth',2)
hold on
semilogx(as,abs(m_div),'x')
xlabel('\alpha')
ylabel('coupling coeff')
legend('curl and cross-wind','div. and along-wind')

f = fit(as',abs(m_curl)','exp1');
st1 = '\alpha';
title(sprintf('exp fit: %2.2f exp ( %2.2f $$%s$$)',f.a,f.b,st1),'interpreter','latex')
set(gca,'fontsize',18)



%% calculate heat fluxes

% thermo/fluid parameters
c_p_air = 1000; % J / kg / K
DT = 0.5;       % K
RH = 0.8;
Lv = 2.26e6;

a = 0;
CD_ref = 1e-3;
CD = (1+a*T_prime)*CD_ref;
QLRH = zeros(N^2,5,2);
QLDT = zeros(N^2,5,2);
for i = 1:2
    if i == 1
        RH_vec = [0.7 0.75 0.8 0.85 0.9];
        RH = [0.7 0.75 0.8 0.85 0.9];
        DT = 0.5;
    elseif i == 2
        RH = 0.8;
        DT_vec = [0.25 0.5 0.75 1 2];
        DT = [0.25 0.5 0.75 1 2];
    end
    for i_RH = 1: length(RH)
        for i_DT = 1:length(DT)
            for j = 1:2
                if j ==1
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
                e_sat = SAM_psatWater(T-DT(i_DT));
                e = RH(i_RH)/100*e_sat;
                r = 0.622 * e ./ max(e, P0-e);
                qa = r./(1+r);
               
                % Q_s(:,:,j) = rho_a.*c_p_air.*CD.*sqrt(U.^2+V.^2).*(DT*ones(size(X)));
                Q_L(:,:,j) = rho_a.*Lv.*CD.*sqrt(U.^2+V.^2).*(qo-qa);
                figure(7)
                subplot(3,2,i+2*(j-1))
                if i == 1
                    QLRH(:,i_RH,j) = reshape(Q_L(:,:,j),N^2,1);
                    %                     plot([1 1]*RH(i_RH),[min(min(Q_L(:,:,j))),max(max(Q_L(:,:,j)))],'o-','linewidth',2)
                    if i_RH == length(RH)
                        boxplot(QLRH(:,:,j),RH(1:i_RH))
                        xlabel('RH')
                    end
                else
                    QLDT(:,i_DT,j) = reshape(Q_L(:,:,j),N^2,1) ;
                    %                     plot([1 1]*DT(i_DT),[min(min(Q_L(:,:,j))),max(max(Q_L(:,:,j)))],'o-','linewidth',2)
                    if i_DT == length(DT)
                        boxplot(QLDT(:,:,j),DT(1:i_DT))
                        xlabel('DT')
                    end
                end
                title([str ' Q_L W m^{-2}'])
                set(gca,'fontsize',18)
            end
        end
    end
    
end



figure(7)
subplot(3,2,5)
boxplot(QLRH(:,:,2)-QLRH(:,:,1),RH_vec)
xlabel('RH')
title(['\Delta Q_L W m^{-2}'])

subplot(3,2,6)
boxplot(QLDT(:,:,2)-QLDT(:,:,1),DT_vec)
xlabel('DT')
title(['\Delta Q_L W m^{-2}'])





