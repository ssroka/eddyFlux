ccc

T1 = 2*pi; % period of signal 1
T2 = pi;   % period of signal 1

L = 500; % length of signal

NT = 3;
N = L;
cT = 2/3; % cutoff period;


%%

t = linspace(0,2*pi*NT,L);

f1 = 2*pi/T1;
f2 = 2*pi/T2;

s1 = cos(f1*t);
s2 = cos(f2*t);

Z = s1+s2;
cf = 1/cT;
% subplot(3,2,1)
% plot(t,s1,'linewidth',2,'displayname',sprintf('signal 1, f = %2.2f rad/s',f1))
% hold on
% plot(t,s2,'linewidth',2,'displayname',sprintf('signal 2, f = %2.2f rad/s',f2))
% plot(t,s,'linewidth',2,'displayname',sprintf('signal'))



Fs = L/NT; % samples per period, L samples for NT periods
F = [-N/2:N/2-1]/N*Fs;



Zhat = fft(Z);







subplot(1,2,1)
plot(t,Z,'linewidth',2);
hold on 
plot(t, cos(2*t),'--','linewidth',2);


P2 = abs(Zhat/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(1,2,2)
% f2 = [fliplr(f) f(2:end)];
% plot(f2,ifftshift(fftshift(P2)))

Zhat = fftshift(abs(Zhat));

plot(F,Zhat,'linewidth',2)

inds = abs(F)<cf;

Zhat(inds) = 0.0;

hold on 

plot(F,Zhat,'g--','linewidth',2)



Z2 = ifft(ifftshift(Zhat));
subplot(1,2,1)
hold on
plot(t,Z2,'o','linewidth',2)


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