
clear;close all;clc
addpath('/Users/ssroka/Documents/MATLAB/util/')

%% ---------------------- User input --------------------------------------
sst_hi = 300; % K
sst_lo = 270; % K
p2p_SST = 8;
uT_coeff = 0.44;
L_km = 350; % km

filter_length_km = 400;  % km 250 km default

% 1 degree of latitude is 100km
% range to investigate is 50km to 400km averaging box
% domain to consider is 
% 4000 km x 4000 km 
xlim_km = [0 10000];

N = 801;
DT = 0.5;
p = 101325; % Pa @ surface
RH = 80; % 
u_bar = 4; % m/s


%% --------------------- Begin --------------------------------------------
xlim = xlim_km*1e3;
L = L_km*1e3;
filter_length = filter_length_km*1e3;

x = linspace(xlim(1),xlim(2),N)';
dx = x(2)-x(1);
N_seg = floor(diff(xlim)/filter_length);
seg_length = floor((N-1)/N_seg); % points/segment
h = floor(seg_length/2); % points/segment

% ---- SST ----
G_window = normpdf(x,mean(xlim),2*L);
SST = G_window./max(G_window).*(0.5*(p2p_SST)*sin(x*2*pi/L));
% SST = (sst_hi-sst_lo)*x/xlim(2)+sst_lo;
plot(SST)
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

F_LAT =max(1e-3*1.2*2.26e6*u.*(qo-qa))
F_SEN =max(1e-3*1.2*4186*u.*(DT*ones(length(qo),1)))

plot(F_LAT)
hold on
plot(F_SEN)
hold off
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


[min(ratio(ratio>0)) max(ratio(ratio>0))]



return

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


%% plot correlations

% 
% figure
% 
% subplot(1,2,1)
% FL_mat = [1.0907e-06 1.3295e-04
% 0.0013     0.0043
% 0.0082     0.0187
% 0.0385     0.0788
% 0.0928     0.1551
%  0.2293     0.2293
% 0.3058     0.5551
% 0.3281     0.6444];
% count = 1;
% for i = 50:50:400
%     plot([i i],FL_mat(count,:),'k-*')
%     count = count + 1;
%     hold on
% end
% 
% xlabel('Filter length [km]')
% title('$$ \frac{\overline{u''(q_o''-q_a'')}}{|\overline{u}|(\overline{q_o}-\overline{q_a})} $$','interpreter','latex')
% 
% set(gcf,'color','w','position',[721     1   720   804])
% set(gca,'fontsize',20)
% % ax1 = gca; % current axes
% % ax1_pos = ax1.Position; % position of first axes
% % ax2 = axes('Position',ax1_pos,...
% %     'XAxisLocation','top',...
% %     'YAxisLocation','right',...
% %     'Color','none');
% 
% count = 1;
% for i = 50:50:400
%     plot([i i]./L_km,FL_mat(count,:),'k-*')
%     count = count + 1;
%     hold on
% end
% 
% 
% subplot(1,2,2)
% u_coeff_rat = [ 0.2   0.0640    0.1017
% 
%  0.25   0.0787    0.1282
% 
%   0.3  0.0928    0.1551
% 
%   0.35  0.1064    0.1824
% 
%   0.4  0.1196    0.2102
% 
%   0.45 0.1323    0.2389 ]
% 
% for i = 1:length(u_coeff_rat)
%     plot([1 1]*u_coeff_rat(i,1),u_coeff_rat(i,2:3),'b-o')
%     hold on
% end
% xlabel('coupling coeff')
% title('$$ \frac{\overline{u''(q_o''-q_a'')}}{|\overline{u}|(\overline{q_o}-\overline{q_a})} $$','interpreter','latex')
% 
% set(gcf,'color','w','position',[721     1   720   804])
% set(gca,'fontsize',20)






