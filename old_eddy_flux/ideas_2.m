% just replace U in the calculation of the ratio with C_D U
clear;close all;clc
for run_ctrl = [1 0]
try
    colormap(othercolor('RdBu9'))
catch
    addpath('~/Documents/MATLAB/util/othercolor/')
    colormap(othercolor('RdBu9'))
end
%% try to repeat Foussard et al. (2019) conclusions with the same parameters


% from table 2

L_km = 9216; % km
y_sst_km = 4500; % km
l_sst_km = 1000; % km
l_an_km = 1500; % km
SST_front = 285; % K
D_SST = 20;% K
D_theta = 10; % K
P0 = 101300; % hPa

eddy_size = [100 500]; % km

NT = 1; % number of time points

%% create domain
y_sst = y_sst_km*1000;
l_sst = l_sst_km*1000;

L = L_km * 1000;
N = 200;
x = linspace(0,L_km,N);
[X, Y] = meshgrid(x,x);


%% create SST fields
l_an = l_an_km*1000;

SST_CTRL = SST_front-D_SST/2*tanh((Y-y_sst_km)/l_sst_km);

SST_eddies_T = zeros(N,N,NT);
if ~run_ctrl
    for j = 1:NT
        for i = 1:300
            x_c= rand*L_km;
            y_c= randn*y_sst_km/3+y_sst_km;
            sig_c = rand*diff(eddy_size)+eddy_size(1);
            SST_eddies_T(:,:,j) = SST_eddies_T(:,:,j) + exp((-(X-x_c).^2-(Y-y_c).^2)/2/sig_c^2);
        end

        % remove zonal mean
        SST_eddies_T(:,:,j) = SST_eddies_T(:,:,j) - repmat(mean(SST_eddies_T(:,:,j),2),1,N);
        a = fzero(@(a) std(reshape(SST_eddies_T(:,:,j),N^2,1)*a)-3,4);
        SST_eddies_T(:,:,j) = SST_eddies_T(:,:,j)*a;
        
        G = exp(-(Y-y_sst_km).^2./l_an_km.^2);
        SST_EDDY_T(:,:,j) = SST_CTRL + SST_eddies_T(:,:,j).*G;
    end
    figure(1)
    contourf(x,x,SST_EDDY_T(:,:,1),30)
    colorbar
    
    SST_EDDY = SST_EDDY_T(:,:,1);
else
    SST_EDDY = SST_CTRL;
end
%% recreate Foussard Fig 3

% a
SST_anomaly = SST_eddies_T(:,:,1);


% d
rho_a = 1.2; % kg m^-3
c_p_air = 1000; % J / kg / K
du = 4;
CH = 1e-3;
DT = 0.5; % K
uT_coeff = 0.44;
u_bar = 4;
u = u_bar/(L_km/2)*(Y-L_km/2)+uT_coeff*SST_anomaly;
RH = 0.8;
Lv = 2.26e6;


%% rec. Fig 4
% ---- qo ----
try
    qo = SAM_qsatWater(SST_EDDY, P0) ;
catch
    addpath('~/Documents/MATLAB/util/')
    qo = SAM_qsatWater(SST_EDDY, P0) ;
end
% ---- qa ----
e_sat = SAM_psatWater(SST_EDDY-DT);
e = RH/100*e_sat;
r = 0.622 * e ./ max(e, P0-e);
qa = r./(1+r);

% F_LAT =max(max(1e-3*1.2*Lv*u.*(qo-qa)))
% F_SEN =max(max(1e-3*1.2*c_p_air*u.*(DT*ones(length(qo),1))))

%% compute fluxes
Q_s = mean(rho_a.*c_p_air.*CH.*abs(u).*(DT*ones(size(X))),2);
Q_L = mean(rho_a.*Lv.*CH.*abs(u).*(qo-qa),2);
figure(2)
if run_ctrl
    Q_s_ctrl = Q_s;
    Q_L_ctrl = Q_L;
    plot(x,Q_s,'k--','linewidth',2,'displayname','CTRL Sensible W/m^2')
    hold on
    plot(x,Q_L,'r--','linewidth',2,'displayname','CTRL Latent W/m^2')
    legend('-dynamiclegend')
    xlabel('km')
else
    Q_s_eddy = Q_s;
    Q_L_eddy = Q_L;
    plot(x,Q_s,'k','linewidth',2,'displayname','EDDY Sensible W/m^2')
    hold on
    plot(x,Q_L,'r','linewidth',2,'displayname','EDDY Latent W/m^2')
    legend('-dynamiclegend')
    xlabel('km')
end
end
ylabel('[W/m^2]')
yyaxis right
plot(x,Q_L_eddy-Q_L_ctrl,'b','displayname','\Delta Q_L')
plot(x,Q_s_eddy-Q_s_ctrl,'b-.','displayname','\Delta Q_L')
ylabel('\Delta Q [W / m^2]')


%% compute stress
CD_ref = 1e-3;
T_prime = SST_anomaly;
alpha = 1;
CD = (1+alpha*T_prime)*CD_ref;
tau = rho_a.*CD.*abs(u.^2);
figure(3)
subplot(2,2,[1 3])
imagesc(tau)
colorbar
title('\tau')
set(gca,'fontsize',18)

% compute the curl
subplot(2,2,2)
c = curl(X,Y,tau,0*tau);
imagesc(c); colorbar
title('\nabla x \tau')
set(gca,'fontsize',18)

% compute the div
subplot(2,2,4)
d = divergence(X,Y,tau,0*tau);
imagesc(d); colorbar
title('\nabla \cdot \tau')
set(gca,'fontsize',18)

%% compute 










