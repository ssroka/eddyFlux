clear;close all;clc

%% free parameters

% air-sea jump for qo
DT = 0.5;% K

p_surf = 101325; % Pa
p_aloft = 95000; % Pa

RH = 0.8;

u_sst_coupling_coeff = 0.44; % accoding to Chelton slides [0.28 0.44]

L = 10; % number of grid cells to make a square over
Lx = 10;
Ly = 10;
h = 0.1;
N = floor(20/h)+1;


%% ---------------- Begin --------------------------------
[sst] = createSST_field(N,Lx,Ly,h);
subplot(1,2,1)
imagesc(sst);colorbar
[sst_bar, sst_prime]  = calc_anomaly(sst,L);

po   = p_surf*ones(size(sst));

[qo] = calc_qo(sst,po,DT);

[qo_bar, qo_prime]  = calc_anomaly(qo,L);

pa   = p_aloft*ones(size(sst));

[qa] = calc_qa(sst,p_surf,p_aloft,DT,RH);

[qa_bar, qa_prime] = calc_anomaly(qa,L);

[u] = calc_u(sst_prime,u_sst_coupling_coeff);

[u_bar, u_prime]   = calc_anomaly(u,L);

numerator = u_prime.*(qo_prime-qa_prime);

[numerator_bar, ~]   = calc_anomaly(numerator,L);

ratio = numerator_bar./(u_bar.*(qo_bar-qa_bar));
subplot(1,2,2)
contourf(ratio);colorbar

[c,h] = contourf(ratio);
clabel(c,h)
set(h,'levellist',[-1:0.01:1])
colorbar

















