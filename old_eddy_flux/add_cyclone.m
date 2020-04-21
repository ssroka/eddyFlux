clear;close all;clc
%% user input

use_uT_coeff_flag = 0;
add_XT_cyc = 1;
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
hold on
try
    colormap(othercolor('RdBu9'))
catch
    addpath('~/Documents/MATLAB/util/othercolor/')
    colormap(othercolor('RdBu9'))
end
%% create wind field
rho_a = 1.2;    % kg m^-3
u_bar = 10;    % m/s

if use_uT_coeff_flag
    uT_coeff = 0.2; % 0.2 - 0.44
    U = u_bar/(L_km/2)*(Y-L_km/2)+uT_coeff*T_prime;
    V = zeros(size(X));
elseif add_XT_cyc
    U = -u_bar*2*(Y-L_km/2)/L_km.*(abs((Y-L_km/2))<1000).*(abs((X-L_km/2))<1000)+...
        u_bar.*((abs((Y-L_km/2))>1000)|(abs((X-L_km/2))>1000));
    V = u_bar*10*(X-L_km/2)/L_km.*(abs((X-L_km/2))<1000).*(abs((Y-L_km/2))<1000);
    P0 = (99700+(1030000-99700)*sqrt((X-L_km/2).^2+(Y-L_km/2).^2)/L_km.*(abs((Y-L_km/2))<1000).*(abs((X-L_km/2))<1000))+...
        1030000.*((abs((Y-L_km/2))>1000)|(abs((X-L_km/2))>1000));
else
    U = u_bar*ones(size(X));
    V = zeros(size(X));
end
% q = quiver(X(1:2:N,1:2:N),Y(1:2:N,1:2:N),U(1:2:N,1:2:N),V(1:2:N,1:2:N));
q = quiver(X,Y,U,V);
set(gca,'fontsize',18)
set(q,'color','k')
figure
contourf(P0)
colorbar
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
    
    [gradSSTX,gradSSTY] = gradient(T_prime,(x(2)-x(1))/100);
    
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
    m = polyfit(gradSSTY(:),cz(:)*1e7,1);
    m_curl(i) = m(1);
    title(['slope = ' num2str(m_curl(i))])
    set(gca,'fontsize',18)
    
    subplot(1,2,2)
    plot(gradSSTX(:),d(:)*1e7,'o')
    xlabel('along-wind SST gradient (deg C / 100 km)')
    ylabel('divergence (N m^-3 x 10^7)')
    m = polyfit(gradSSTX(:),d(:)*1e7,1);
    m_div(i) = m(1);
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
    qo = SAM_qsatWater(T, P0*.9) ;
catch
    addpath('~/Documents/MATLAB/util/')
    qo = SAM_qsatWater(T, P0*.9) ;
end
% ---- qa ----
e_sat = SAM_psatWater(T-DT);
e = RH/100*e_sat;
r = 0.622 * e ./ max(e, P0-e);
qa = r./(1+r);

a = 0.25;
CD_ref = 1e-3;
CD = (1+a*T_prime)*CD_ref;

Q_s(:,:,j) = rho_a.*c_p_air.*CD.*sqrt(U.^2+V.^2).*(DT*ones(size(X)));
Q_L(:,:,j) = rho_a.*Lv.*CD.*sqrt(U.^2+V.^2).*(qo-qa);

figure(7)
subplot(3,2,1+2*(j-1))
imagesc(x,x,Q_s(:,:,j))
colorbar
xlabel('km')
ylabel('km')
title([str ' Q_s W m^{-2}'])
set(gca,'fontsize',18)

subplot(3,2,2+2*(j-1))
imagesc(x,x,Q_L(:,:,j))
colorbar
xlabel('km')
ylabel('km')
title([str ' Q_L W m^{-2}'])
set(gca,'fontsize',18)
end

subplot(3,2,5)
imagesc(x,x,Q_s(:,:,2)-Q_s(:,:,1))
colorbar
xlabel('km')
ylabel('km')
title(['\Delta Q_s W m^{-2}'])
set(gca,'fontsize',18)
subplot(3,2,6)
imagesc(x,x,Q_L(:,:,2)-Q_L(:,:,1))
colorbar
xlabel('km')
ylabel('km')
title(['\Delta Q_L W m^{-2}'])
set(gca,'fontsize',18)
