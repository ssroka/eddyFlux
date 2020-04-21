
clear;close all;clc

DT = 0.5;
P0  = 101300;    % Pa (surface pressure)
Ts = 280:312;
for i =1:length(Ts)
    T = Ts(i);
    try
        qo(i) = SAM_qsatWater(T, P0) ;
    catch
        addpath('~/Documents/MATLAB/util/')
        qo(i) = SAM_qsatWater(T, P0) ;
    end
end
plot(Ts,qo,'o')
hold on
p = polyfit(Ts,qo,4);
plot(Ts,polyval(p,Ts),'r')
figure
T = Ts;
dq_star_dT = 3*T.^2*p(1)+2*T.^1*p(2)+p(3);
plot(Ts,dq_star_dT)
%
% % ---- qa ----
% e_sat = SAM_psatWater(T-DT);
% e = RH/100*e_sat;
% r = 0.622 * e ./ max(e, P0-e);
% qa = r./(1+r);