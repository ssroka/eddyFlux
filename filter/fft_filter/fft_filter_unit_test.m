% =========================================================================
% fft_filter_unit_test
% =========================================================================
% Created By    : Sydney Sroka
% Last Edited By: Sydney Sroka
% Date          : 07/14/2020
%
% ------------------------
% This script was written to practice using fft to filter a 1D signal

clear;close all;clc

% if you change the period length to anything that isn't a multiply of pi
% (or you could change the length of the input NT to some
% non-whole number) the output is not exact

T1 = 3*2*pi; % [s] period of signal 1
T2 = 0.5*2*pi;   % [s] period of signal 1

N = 601; % length of signal

NT = 3.5;

cT = 2*2*pi; % cutoff period;

%%

t = linspace(0,2*pi*NT,N);

f1 = 1/T1;
f2 = 1/T2;

cf = 1/cT;

s1 = cos(2*pi/T1*t);
s2 = cos(2*pi/T2*t);

s = s1.*s2;

subplot(4,2,1)
plot(t,s1,'linewidth',2,'displayname',sprintf('signal 1, f = %2.2f rad/s',f1))
hold on
plot(t,s2,'linewidth',2,'displayname',sprintf('signal 2, f = %2.2f rad/s',f2))
plot(t,s,'linewidth',2,'displayname',sprintf('signal'))


Fs = N/(2*pi*NT); % samples per period, L samples for NT periods
F = [-N/2:N/2-1]/N*Fs;
tr = s-detrend(s);
s = detrend(s);
Zhat = fft(s);

Zhat = fftshift(abs(Zhat));

Zhat_LP = Zhat;
Zhat_HP = Zhat;

inds = abs(F)>cf;
Zhat_LP(inds) = 0.0;

inds = abs(F)<cf;
Zhat_HP(inds) = 0.0;

Z_LP = ifft(ifftshift(Zhat_LP))+tr;
Z_HP = ifft(ifftshift(Zhat_HP));


subplot(4,2,2)
plot(F,abs(Zhat),'linewidth',2,'displayname','FFT(s)');
hold on 

subplot(4,2,3)
plot(F,Zhat,'linewidth',2,'displayname','FFT(s)');
hold on 
plot(F,abs(Zhat_LP),'--','linewidth',2,'displayname','FFT(s)');
title('low pass')

subplot(4,2,4)
plot(F,Zhat,'linewidth',2,'displayname','FFT(s)');
hold on 
plot(F,abs(Zhat_HP),'--','linewidth',2,'displayname','FFT(s)');
title('high pass')


plt_ind = 1:5:length(s1);

subplot(4,2,5)
plot(t,s1,'linewidth',2,'displayname',sprintf('signal 1, f = %2.2f rad/s',f1))
hold on
plot(t,s2,'linewidth',2,'displayname',sprintf('signal 2, f = %2.2f rad/s',f2))
plot(t,s,'linewidth',2,'displayname',sprintf('signal'))
plot(t(plt_ind),Z_LP(plt_ind),'o','linewidth',2,'displayname','FFT(s)');
title('low pass')

subplot(4,2,6)
plot(t,s1,'linewidth',2,'displayname',sprintf('signal 1, f = %2.2f rad/s',f1))
hold on
plot(t,s2,'linewidth',2,'displayname',sprintf('signal 2, f = %2.2f rad/s',f2))
plot(t,s,'linewidth',2,'displayname',sprintf('signal'))
plot(t(plt_ind),Z_HP(plt_ind),'o','linewidth',2,'displayname','FFT(s)');
title('high pass')

plot(t(plt_ind),Z_HP(plt_ind)+Z_LP(plt_ind),'go','linewidth',2,'displayname','FFT(s)');

subplot(4,2,[7 8])
plot(t(plt_ind),Z_LP(plt_ind)+Z_HP(plt_ind)-s(plt_ind))

set(gcf,'position',[440           1        1001         797])


% s = [zeros(1,50) s(51:450) zeros(1,50)];

% [X, XHP, f, y, y2] = fftf(t, s, cf);
% close all
% figure
% subplot(4,1,1)
% plot(t,s)
% subplot(4,1,2)
% plot(t,X)
% hold on
% plot(t,s1)
% subplot(4,1,3)
% plot(t,XHP)
% hold on
% plot(t,s2)
% subplot(4,1,4)
% plot(t,s)
% hold on
% plot(t,X+XHP)
% plot(t,s-(X+XHP),'linewidth',3)

% plot(f,y(1:1+length(f)/2))
% plot(f,y2(1:1+length(f)/2))

% subplot(1,2,2)
% f2 = [fliplr(f) f(2:end)];
% plot(f2,ifftshift(fftshift(P2)))



% subplot(1,2,1)
% hold on
% plot(t,Z2,'--','linewidth',2)




% subplot(1,2,2)
% 
% Zhat = fft2(Z);
% % P2 = abs(Zhat./L);
% % P1 = P2(1:L/2+1,1:L/2+1);
% % P1(2:end-1,2:end-1) = 2*P1(2:end-1,2:end-1);
% % f = Fs*(0:(L/2))/L;
% contourf(abs(fftshift(Zhat)));
% set(gca,'ydir','normal')

% [X,Y] = meshgrid(x,x);
% 
% P = peaks(20);
% Z = repmat(P,[5 10]);

% Z = cos(X)+cos(Y-pi/4)+rand(L,L);
% P2 = abs(Zhat/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);