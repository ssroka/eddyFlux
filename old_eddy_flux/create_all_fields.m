clear;close all;clc
%% User Input

x_range = [-10 10];
y_range = [-10 10];

h = 0.1;

% characteristic eddy length scale
Lx       = 4;
Ly       = 4;

sst_hi = 300;
sst_lo = 290;

u_hi = 7;
u_lo = 2;

q_o_hi = 0.030; %kg/kg
q_o_lo = 0.020; %kg/kg

q_a_hi = 0.015; %kg/kg
q_a_lo = 0.005;  %kg/kg

%% Begin

% ============== set SST ==========================================

N = floor(diff(x_range)/h)+1;
dx = round(Lx/h);
dy = round(Ly/h);
% create eddy field
sst = zeros(N);

mu_x = 0;
mu_y = 0;
sig_x = 0.5;
sig_y = 0.7;

x = repmat(linspace(-Lx/2,Lx/2,dx),dy,1);
y = repmat([linspace(-Ly/2,Ly/2,dy)]',1,dx);
alpha = 0.1;
filter_flag = true;

[sst_mag,sstf_mag]=create_field(N,dx,dy,sig_x,sig_y,mu_x,mu_y,x,y,filter_flag,alpha);
sst = sst_mag*(sst_hi-sst_lo)+sst_lo; % 
sstf = sstf_mag*(sst_hi-sst_lo)+sst_lo; % 

% subplot(2,2,1)
% imagesc(sst)
% colorbar
% subplot(2,2,2)
% imagesc(sstf)
% colorbar
% subplot(2,2,3)
% imagesc(sst-sstf)
% colorbar
% ============== set q_o ==========================================

% q should be somewhat spatialy correlated with SST

% q_o

% u should be somewhat spatialy correlated with SST
q_o_mu_x = 0;
q_o_mu_y = 0;
q_o_sig_x = 0.25;
q_o_sig_y = 0.35;
alpha = 0.1;

[q_o_mag,q_of_mag]=create_field(N,dx,dy,q_o_sig_x,q_o_sig_y,q_o_mu_x,q_o_mu_y,x,y,filter_flag,alpha);

q_o = q_o_mag*(q_o_hi-q_o_lo)+q_o_lo; % 
q_of = q_of_mag*(q_o_hi-q_o_lo)+q_o_lo; % 

% subplot(2,2,4)
% imagesc(q_o)
% colorbar
subplot(2,2,1)
imagesc(q_of)
title('$$q_o$$ [kg/kg]','interpreter','latex')
colorbar
% subplot(2,2,6)
% imagesc(q_o-q_of)
% colorbar

% ============== set q_a ==========================================

% q should be somewhat spatialy correlated with SST

% q_a

% u should be somewhat spatialy correlated with SST
q_a_mu_x = 0;
q_a_mu_y = 0;
q_a_sig_x = 0.25;
q_a_sig_y = 0.35;
alpha = 0.1;

[q_a_mag,q_af_mag]=create_field(N,dx,dy,q_a_sig_x,q_a_sig_y,q_a_mu_x,q_a_mu_y,x,y,filter_flag,alpha);

q_a = (q_a_mag-min(q_a_mag(:)))*(q_a_hi-q_a_lo)+q_a_lo; % 
q_af = (q_af_mag-min(q_af_mag(:)))*(q_a_hi-q_a_lo)+q_a_lo; % 

% subplot(2,2,7)
% imagesc(q_a)
% colorbar
subplot(2,2,2)
imagesc(q_af)
title('$$q_a$$ [kg/kg]','interpreter','latex')
colorbar
% subplot(2,2,9)
% imagesc(q_a-q_af)
% colorbar



% ============== set u ==========================================

% u should be somewhat spatialy correlated with SST
u_mu_x = 0;
u_mu_y = 0;
u_sig_x = 0.2;
u_sig_y = 0.3;
alpha = 0.1;

[u_mag,uf_mag]=create_field(N,dx,dy,u_sig_x,u_sig_y,u_mu_x,u_mu_y,x,y,filter_flag,alpha);
u = u_mag*(u_hi-u_lo)+u_lo; % 
uf = uf_mag*(u_hi-u_lo)+u_lo; % 

% subplot(2,2,10)
% imagesc(u)
% colorbar
subplot(2,2,3)
imagesc(u)
title('$$u [m/s] $$','interpreter','latex')
colorbar
% subplot(2,2,12)
% imagesc(u-uf)
% colorbar


function [e,ef]=create_field(N,dx,dy,sig_x,sig_y,mu_x,mu_y,x,y,filter_flag,alpha)

e = zeros(N);

for i = 1:floor(N/dx)
    X = 1+(i-1)*dx:dx+(i-1)*dx;
    for j = 1:floor(N/dy)
        Y = 1+(j-1)*dy:dy+(j-1)*dy;
        sx = sig_x+rand*0.5;
        sy = sig_y+rand*0.5;
        mx = mu_x+rand*0.5;
        my = mu_y+rand*0.5;
        e(Y,X) = (-1).^(mod(i,2))*(-1).^(mod(j,2))*exp(-(x - mx).^2/2/sx^2-(y - my).^2/2/sy^2);
    end
end

if filter_flag
    % low pass filter
    ef = zeros(size(e));
    ef(:,1) = e(:,1);
    %  alpha = 0.1;
    for i = 2:size(e,1)
        for j = 2:size(e,2)
            ef(i,j) = 0.5*(ef(i-1,j) + alpha*(e(i,j)-ef(i-1,j)) +...
                ef(i,j-1) + alpha*(e(i,j)-ef(i,j-1)) );
        end
    end
end

end












