
clear;close all;clc
addpath('/Users/ssroka/Documents/MATLAB/util/')

%% ---------------------- User input --------------------------------------
sst_hi = 300; % K
sst_lo = 295; % K
uT_coeff = 0.3;
Lx_km = 1000; % km
Ly_km = 8000; % km

filter_length_km = 250;  % km

% 1 degree of latitude is 100km
% range to investigate is 50km to 400km averaging box
% domain to consider is 
% 4000 km x 4000 km 
xlim_km = [0 4000];

N = 401;
DT = 0.5;
p = 101325; % Pa @ surface
RH = 80; % 
u_bar = 4; % m/s


%% --------------------- Begin --------------------------------------------
xlim = xlim_km*1e3;
Lx = Lx_km*1e3;
Ly = Ly_km*1e3;
filter_length = filter_length_km*1e3;

x = linspace(xlim(1),xlim(2),N)';
[X,Y] = meshgrid(x,x);

dx = x(2)-x(1);
N_seg = floor(diff(xlim)/filter_length);
seg_length = floor((N-1)/N_seg); % points/segment
h = floor(seg_length/2); % points/segment

% ---- SST ----
SST = (sst_hi-sst_lo)*0.5*(sin(X*2*pi/Lx).*sin(Y*2*pi/Ly))+sst_lo; %

[dSST_dx,dSST_dy] = gradient(SST,dx); 

quiver(X,Y,dSST_dx,dSST_dy)


SST_prime = zeros(N);

for ipx = h+1:N-h-1
    inds_x = (ipx-h:ipx+h);
    for ipy = h+1:N-h-1
        inds_y = (ipy-h:ipy+h);
        SST_prime(ipy,ipx) = SST(ipy,ipx)-mean(mean(SST(inds_y,inds_x)));
    end
end
% ---- u ----
u = u_bar+SST_prime*uT_coeff;
u_prime = zeros(N);

% ---- qo ----
qo = SAM_qsatWater(SST-DT, p) ; 
qo_prime = zeros(N);

% ---- qa ----
e_sat = SAM_psatWater(SST-DT);
e = RH/100*e_sat;
r = 0.622 * e ./ max(e, p-e);
qa = r./(1+r);
qa_prime = zeros(N);

inds_minh = h+1:N-h-1;
inds_min2h = 2*h+1:N-2*h-1;


for ipx = h+1:N-h-1
    inds_x = (ipx-h:ipx+h);
    for ipy = h+1:N-h-1
        inds_y = (ipy-h:ipy+h);
    SST_prime(ipy,ipx) = SST(ipx)-mean(mean(SST(inds_y,inds_x)));
    u_prime(ipy,ipx)   = u(ipx)-mean(mean(u(inds_y,inds_x)));
    qo_prime(ipy,ipx)  = qo(ipx)-mean(mean(qo(inds_y,inds_x)));
    qa_prime(ipy,ipx)  = qa(ipx)-mean(mean(qa(inds_y,inds_x)));
    end
end

ratio = zeros(N);

for ipx = h+1:N-h-1
    inds_x = (ipx-h:ipx+h);
    for ipy = h+1:N-h-1
        inds_y = (ipy-h:ipy+h);
    ratio(ipy,ipx) = mean(mean(u_prime(inds_y,inds_x).*(qo_prime(inds_y,inds_x)-qa_prime(inds_y,inds_x))))/(mean(mean(u(inds_x,inds_y)))*(mean(mean(qo(inds_x,inds_y)))-mean(mean(qa(inds_x,inds_y)))));
    end
end

subplot(3,2,1)
imagesc(x(inds_minh)*1e-3,x(inds_minh)*1e-3,SST_prime(inds_minh,inds_minh))
title('SST'' [K]')
xlabel(' [km] ')
ylabel(' [km] ')
colorbar
set(gca,'fontsize',20)

subplot(3,2,2)
imagesc(x(inds_minh)*1e-3,x(inds_minh)*1e-3,u(inds_minh,inds_minh))
title('u [m/s]')
xlabel(' [km] ')
ylabel(' [km] ')
colorbar
set(gca,'fontsize',20)

addpath('/Users/ssroka/Documents/MATLAB/util/')

subplot(3,2,3)
imagesc(x(inds_minh)*1e-3,x(inds_minh)*1e-3,qo_prime(inds_minh,inds_minh))
xlabel(' [km] ')
ylabel(' [km] ')
title('q_o [kg/kg]')
colorbar
clim_o = get(gca,'clim');
set(gca,'fontsize',20)
subplot(3,2,4)
imagesc(x(inds_minh)*1e-3,x(inds_minh)*1e-3,qa_prime(inds_minh,inds_minh))
xlabel(' [km] ')
ylabel(' [km] ')
title('q_a [kg/kg]')
colorbar
set(gca,'clim',clim_o);

set(gca,'fontsize',20)

subplot(3,2,5:6)
imagesc(x(inds_min2h)*1e-3,x(inds_min2h)*1e-3,ratio(inds_min2h,inds_min2h))
xlabel(' [km] ')
colorbar
title('$$ \frac{\overline{u''(q_o''-q_a'')}}{|\overline{u}|(\overline{q_o}-\overline{q_a})} $$','interpreter','latex')


set(gcf,'color','w','position',[721     1   720   804])
set(gca,'fontsize',20)


[min(ratio(ratio>0)) max(ratio(ratio>0))]



% filter length = 250 km
% u_coeff = 0.3 default
% sst_hi = 303; % K
% sst_lo = 293; % K



% filter length = 50 km  1.0907e-06 1.3295e-04
% filter length = 100 km 0.0013     0.0043
% filter length = 150 km 0.0082     0.0187
% filter length = 200 km 0.0385     0.0788
% filter length = 250 km 0.0928     0.1551
% filter length = 300 km 0.2293     0.2293
% filter length = 350 km 0.3058     0.5551
% filter length = 400 km 0.3281     0.6444



% u_coeff

 % 0.2   0.0640    0.1017

 % 0.25   0.0787    0.1282

  % 0.3  0.0928    0.1551

  % 0.35  0.1064    0.1824

  % 0.4  0.1196    0.2102

  % 0.45 0.1323    0.2389




















