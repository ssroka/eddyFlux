% =========================================================================
% fft_filter_unit_test_2D
% =========================================================================
% Created By    : Sydney Sroka
% Last Edited By: Sydney Sroka
% Date          : 07/14/2020
%
% ------------------------
% This script was written to practice using fft2 to filter a 2D signal

clear;close all;clc

% if you change the period length to anything that isn't a multiply of pi
% (or you could change the length of the input NT to some
% non-whole number) the output is not exact

f1 = 7; % [1/s] frq of signal 1
f2 = 1;  % [1/s] period of signal 2

L = 800; % length of signal
N = L;

NT = 1;

cf = 5; % cutoff frq;

debug_flag = false;
%%

t = linspace(0,2*pi*NT,L);
dt = t(2)-t(1);

[X,Y] = meshgrid(t,t);

s1 = cos(2*pi*f1*X);%+cos(f1*Y);
s2 = cos(2*pi*f2*X);%+cos(f2*Y);

s = s1+s2;

% subplot(4,2,1)
% contourf(t,t,s1)
% title(sprintf('signal 1, f = %2.2f rad/s',f1))
%
% subplot(4,2,2)
% contourf(t,t,s2)
% title(sprintf('signal 2, f = %2.2f rad/s',f2))

% subplot(4,2,3)
% contourf(t,t,s1,'b','displayname',sprintf('signal 1, f = %2.2f rad/s',f1))
% hold on
% contourf(t,t,s2,'r','displayname',sprintf('signal 2, f = %2.2f rad/s',f2))
% contourf(t,t,s,'k','displayname',sprintf('signal'))
% title('signal1 + signal 2')



% [S_LP,S_HP] = FFT2D_filter(s,dt,cf,debug_flag,t,t);
close all

sum_s1s2 = s1+s2;
[sum_LP,sum_HP] = FFT2D_filter(sum_s1s2,dt,cf,debug_flag,t,t);


[sum_HP_LP,sum_HP_HP] = FFT2D_filter(sum_LP,dt,cf,debug_flag,t,t);
[sum_HP_LP,sum_HP_HP] = FFT2D_filter(sum_HP,dt,cf,debug_flag,t,t);

% debug_flag = true;
prod_s1s2 = s1.*s2;
[prod_LP,prod_HP] = FFT2D_filter(prod_s1s2,dt,cf,debug_flag,t,t);

%%
figure(3)
subplot(4,3,1)
contourf(t,t,s1)
title(sprintf('signal 1, f = %2.2f rad/s',f1))

subplot(4,3,4)
contourf(t,t,s2)
title(sprintf('signal 2, f = %2.2f rad/s',f2))

subplot(4,3,7)
plot(t,s1(1,:))
hold on
plot(t,s2(1,:))
title(sprintf('signal 1, f = %2.2f rad/s\\signal 2, f = %2.2f rad/s',f1,f2))

%
subplot(4,3,2)
contourf(t,t,sum_s1s2)
title(sprintf('signal 1 + signal 2'))

subplot(4,3,8)
contourf(t,t,sum_LP)
title(sprintf('LP(signal 1 + signal 2)'))

subplot(4,3,11)
contourf(t,t,sum_HP)
title(sprintf('HP(signal 1 + signal 2)'))

%
subplot(4,3,3)
contourf(t,t,prod_s1s2)
title(sprintf('signal 1 x signal 2'))

subplot(4,3,6)
contourf(t,t,prod_LP)
title(sprintf('LP(signal 1 x signal 2)'))

subplot(4,3,9)
[~,h] = contourf(t,t,prod_HP);
set(h,'edgecolor','none')
title(sprintf('HP(signal 1 x signal 2)'))

subplot(4,3,12)
[~,h] = contourf(t,t,prod_s1s2 -(prod_LP+prod_HP));
set(h,'edgecolor','none')
title(sprintf('s1s2-(LP(s1 x s2)+HP(s1 x s2))\\should vanish'))

for i = 1:12
    subplot(4,3,i)
    set(gca,'xtick',[0:6]*pi,'ytick',[0:6]*pi,'xlim',[0 2*pi])
    colorbar
end
    set(gcf,'position',[132           7        1120         798])
