
clear;close all;clc
addpath('/Users/ssroka/Documents/MATLAB/util/')

%% ---------------------- User input --------------------------------------
save_file_str = 'change_abs_ocean_temp';
sst_lo_vec = 290:305; % K
% sst_hi = 300; % K
uT_coeff = 0.3;
L_km = 600; % km

filter_length_km = 200;  % km 250 km default

% 1 degree of latitude is 100km
% range to investigate is 50km to 400km averaging box
% domain to consider is
% 4000 km x 4000 km
xlim_km = [0 4000];

N = 801;
DT = 0.5;
p = 101325; % Pa @ surface
RH = 80; % 80 
u_bar = 4; % m/s

plot_flag = true;

param = sst_lo_vec;
param_str = '$$T_{ocean}$$ [K]';
%% --------------------- Begin --------------------------------------------
change_var = zeros(length(param),3);

xlim = xlim_km*1e3;
L = L_km*1e3;


for iv = 1:length(param)
sst_lo = sst_lo_vec(iv);
sst_hi = sst_lo+5;
% --filt length
filter_length = filter_length_km*1e3;

x = linspace(xlim(1),xlim(2),N)';
dx = x(2)-x(1);
N_seg = floor(diff(xlim)/filter_length);
seg_length = floor((N-1)/N_seg); % points/segment
h = floor(seg_length/2); % points/segment

% ---- SST ----
SST = 0.5*(sst_hi-sst_lo)*sin(x*2*pi/L)+sst_lo+0.5*(sst_hi-sst_lo);
SST_prime = zeros(N,1);

for ip = h+1:N-h-1
    SST_prime(ip) = SST(ip)-mean(SST(ip-h:ip+h));
end
% ---- u ----
u = u_bar+SST_prime*uT_coeff;
u_prime = zeros(N,1);

% ---- qo ----
qo = SAM_qsatWater(SST, p) ;
qo_prime = zeros(N,1);

% ---- qa ----
e_sat = SAM_psatWater(SST-DT);
e = RH/100*e_sat;
r = 0.622 * e ./ max(e, p-e);
qa = r./(1+r);
qa_prime = zeros(N,1);

inds_minh = h+1:N-h-1;
inds_min2h = 2*h+1:N-2*h-1;

for ip = h+1:N-h-1
    inds = (ip-h:ip+h);
    SST_prime(ip) = SST(ip)-mean(SST(inds));
    u_prime(ip)   = u(ip)-mean(u(inds));
    qo_prime(ip)  = qo(ip)-mean(qo(inds));
    qa_prime(ip)  = qa(ip)-mean(qa(inds));
end

for ip = 2*h+1:N-2*h-1
    inds = (ip-h:ip+h);
    ratio(ip) = mean(u_prime(inds).*(qo_prime(inds)-qa_prime(inds)))/(mean(u(inds))*(mean(qo(inds))-mean(qa(inds))));
end

change_var(iv,:) = [(param(iv)) min(ratio(ratio>0)) max(ratio(ratio>0))];

end

save(save_file_str,'change_var','param_str')





% plot
if plot_flag
    
    
    subplot(3,1,1)
    yyaxis left
    plot(x(inds_minh)*1e-3,SST(inds_minh),'linewidth',2)
    ylabel('SST'' [K]')
    xlabel(' [km] ')
    hold on
    yyaxis right
    plot(x(inds_minh)*1e-3,u(inds_minh),'-','linewidth',2)
    ylabel('u [m/s]')
    set(gca,'fontsize',20)
    
    addpath('/Users/ssroka/Documents/MATLAB/util/')
    
    subplot(3,1,2)
    plot(x(inds_minh)*1e-3,qo_prime(inds_minh),'linewidth',2)
    hold on
    plot(x(inds_minh)*1e-3,qa_prime(inds_minh),'linewidth',2)
    xlabel(' [km] ')
    legend('q_o'' [kg/kg]','q_a'' [kg/kg]')
    set(gca,'fontsize',20)
    
    subplot(3,1,3)
    plot(x(inds_min2h)*1e-3,ratio(inds_min2h),'linewidth',2)
    xlabel(' [km] ')
    title('$$ \frac{\overline{u''(q_o''-q_a'')}}{|\overline{u}|(\overline{q_o}-\overline{q_a})} $$','interpreter','latex')
    
    
    set(gcf,'color','w','position',[721     1   720   804])
    set(gca,'fontsize',20)
    
    
    
end
