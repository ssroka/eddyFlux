ccc

Tfinal = 2*pi*2;
N = 1000;
t = linspace(0,Tfinal,N);
x = cos(t)+cos(8*t);
x(500:520) = NaN;

for i = 1:8
subplot(4,2,i)
plot(t,x)
hold on
end

%% test low pass filter

subplot(4,2,1)
title('low pass')

% cos(t) is 1 cycle every 2*pi which is approx 1/6
% cos(8 t) is 8 cycles every 2*pi which is approx 8/6

% low pass fileter out the fast cycle by setting the cutoff below the fast
% cycle and above the slow cycle at 4/6
dt = t(2)-t(1);
cf = 4/6;
xf = lanczosfilter(x,dt,cf);
subplot(4,2,3)
plot(t,xf,'--','linewidth',2)

% low pass fileter out the slow and fast cycles by setting the cutoff below
% the slow cycle at 0.5/6
dt = t(2)-t(1);
cf = 0.5/6;
xf = lanczosfilter(x,dt,cf);
subplot(4,2,5)
plot(t,xf,'--','linewidth',2)

% t is 100 points, so the highest frequency cycle that could be detected 
% is a cycle that happens 50 time between 0 and Tfinal - multiply by 0.99
% to avoid equalling the Nyquist
dt = t(2)-t(1);
cf = (N/2)/Tfinal*.99;
xf = lanczosfilter(x,dt,cf);
subplot(4,2,7)
plot(t,xf,'--','linewidth',2)


%% test high pass filter

% cos(t) is 1 cycle every 2*pi which is approx 1/6
% cos(8 t) is 8 cycles every 2*pi which is approx 8/6

subplot(4,2,2)
title('high pass')


% low pass fileter out the fast cycle by setting the cutoff below the fast
% cycle and above the slow cycle at 4/6
dt = t(2)-t(1);
cf = 4/6;
xf = lanczosfilter(x,dt,cf,[],'high');
subplot(4,2,4)
plot(t,xf,'--','linewidth',2)

% low pass fileter out the slow and fast cycles by setting the cutoff below
% the slow cycle at 0.5/6
dt = t(2)-t(1);
cf = 0.5/6;
xf = lanczosfilter(x,dt,cf,[],'high');
subplot(4,2,6)
plot(t,xf,'--','linewidth',2)

% t is 100 points, so the highest frequency cycle that could be detected 
% is a cycle that happens 50 time between 0 and Tfinal - multiply by 0.99
% to avoid equalling the Nyquist
dt = t(2)-t(1);
cf = (N/2)/Tfinal*.99;
xf = lanczosfilter(x,dt,cf,[],'high');
subplot(4,2,8)
plot(t,xf,'--','linewidth',2)



set(gcf,'position',[25           1        1416         787])







