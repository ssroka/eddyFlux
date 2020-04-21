clear;close all;clc
%{
Default
sst_hi = 300; % K
sst_lo = 295; % K
uT_coeff = 0.3;
L_km = 600; % km
filter_length_km = 200;  % km
RH = 80
DT = 0.5;
u_bar = 4; % m/s
%}

for m = 1:7
    subplot(2,4,m)
 switch m
    case 1
        load change_RH
    case 2
        load change_uT
    case 3
        load change_DSST
    case 4
        load change_filt_len
    case 5
        load change_DT
    case 6
        load change_ubar
     case 7
         load change_abs_ocean_temp
end

for i = 1:size(change_var,1)
    plot([1 1]*change_var(i,1),change_var(i,2:3),'b-*')
    hold on
end
xlabel(param_str,'interpreter','latex')
title('$$ \frac{\overline{u''(q_o''-q_a'')}}{|\overline{u}|(\overline{q_o}-\overline{q_a})} $$','interpreter','latex')

set(gcf,'color','w','position',[721     1   720   804])
set(gca,'fontsize',20)
end