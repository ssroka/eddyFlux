
clear;close all;clc
createSST_field
u_prime = u-mean(u(:));
q_o_prime = q_o - mean(q_o(:));
q_a_prime = q_a - mean(q_a(:));

rho = 1.2;
CD = 1e-3;

filtered_result_num   = u_prime.*(q_o_prime-q_a_prime);
filtered_result_denom = abs(mean(u(:))).*(mean(q_o(:))-mean(q_a(:)));


subplot(2,2,4)
imagesc(filtered_result_num/filtered_result_denom)
title('$$ \frac{\overline{u''(q_o''-q_a'')}}{|\overline{u}|(\overline{q_o}-\overline{q_a})} $$','interpreter','latex')
colorbar




