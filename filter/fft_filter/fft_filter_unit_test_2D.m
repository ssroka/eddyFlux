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

T1 = 2*pi; % [s] period of signal 1
T2 = pi;   % [s] period of signal 1

L = 500; % length of signal
N = L;

NT = 3;

cT = 1.5*pi; % cutoff period;

%%

t = linspace(0,2*pi*NT,L);

f1 = 1/T1;
f2 = 2*pi/T2;

cf = 1/cT;

[X,Y] = meshgrid(t,t);

s1 = cos(2*pi/T1*X)+cos(2*pi/T1*Y);
s2 = cos(2*pi/T2*X)+cos(2*pi/T2*Y);

s = s1+s2;

subplot(4,2,1)
contourf(t,t,s1)
title(sprintf('signal 1, f = %2.2f rad/s',f1))

subplot(4,2,2)
contourf(t,t,s2)
title(sprintf('signal 2, f = %2.2f rad/s',f1))

subplot(4,2,3)
contourf(t,t,s1,'b','displayname',sprintf('signal 1, f = %2.2f rad/s',f1))
hold on
contourf(t,t,s2,'r','displayname',sprintf('signal 2, f = %2.2f rad/s',f2))
contourf(t,t,s,'k','displayname',sprintf('signal'))
title('signal1 + signal 2')


Fs = L/(2*pi*NT); % samples per period, L samples for NT periods
F = [-N/2:N/2-1]/N*Fs;

[Xhat,Yhat] = meshgrid(F,F);



Zhat = fft2(s);

Zhat = fftshift(abs(Zhat));

Zhat_LP = Zhat;
Zhat_HP = Zhat;

inds = sqrt(Xhat.^2+Yhat.^2)>cf;
Zhat_LP(inds) = 0.01;

inds = sqrt(Xhat.^2+Yhat.^2)<cf;
Zhat_HP(inds) = 0.01;

Z_LP = ifft2(ifftshift(Zhat_LP));
Z_HP = ifft2(ifftshift(Zhat_HP));


subplot(4,2,4)
contourf(F,F,abs(Zhat),'displayname','FFT(s)');
hold on 
title('FFT(signal1 + signal 2)')


subplot(4,2,5)
contourf(F,F,abs(Zhat_LP),'--','displayname','FFT(s)');
title('low pass')

subplot(4,2,6)
contourf(F,F,abs(Zhat_HP),'--','displayname','FFT(s)');
title('high pass')

subplot(4,2,7)
contourf(t,t,abs(Z_LP))
title('low pass')

subplot(4,2,8)
contourf(t,t,abs(Z_HP))
title('high pass')

set(gcf,'position',[440           1        1001         797])

