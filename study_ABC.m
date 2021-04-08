














CD = rho_a*0.0014;

ax = subplot(2,4,1);
[~,h] = contourf(lon_plot,lat_plot,CD*A(:,:,count)');
title('A:$\rho_a C_D \overline{U}\overline{\Delta h}$','interpreter','latex')
format_fig(h,ax)

 ax = subplot(2,4,2);
    [~,h] = contourf(lon_plot,lat_plot,CD*B(:,:,count)');
    title('B:$\rho_a C_D\overline{U''\Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,3);
    [~,h] = contourf(lon_plot,lat_plot,CD*C1(:,:,count)');
    title('C1:$\rho_a C_D\overline{T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,4);
    [~,h] = contourf(lon_plot,lat_plot,CD*C2(:,:,count)');
    title('C2:$\rho_a C_D\overline{U''T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,5);
    [~,h] = contourf(lon_plot,lat_plot,CD*C3(:,:,count)');
    title('C3:$\rho_a C_D\overline{U'' T_o''\alpha \Delta h''}$','interpreter','latex')
    format_fig(h,ax)
    
     ax = subplot(2,4,6);
    [~,h] = contourf(lon_plot,lat_plot,CD*D(:,:,count)');
    title('D:$\rho_a C_D\overline{U''\overline{\Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,7);
    [~,h] = contourf(lon_plot,lat_plot,CD*E1(:,:,count)');
    title('E1:$\rho_a C_D\overline{\overline{U} T_o''\overline{\alpha \Delta h}}$','interpreter','latex')
    format_fig(h,ax)
    
    ax = subplot(2,4,8);
    [~,h] = contourf(lon_plot,lat_plot,CD*E2(:,:,count)');
    title('E2:$\rho_a C_D\overline{\overline{U}  \Delta h''}$','interpreter','latex')
    format_fig(h,ax)


    figure
    ax = subplot(1,4,1);
    [~,h] = contourf(lon_plot,lat_plot,h_diff');
    title('$\Delta h$','interpreter','latex')
    format_fig(h,ax)
    ax = subplot(1,4,2);
    [~,h] = contourf(lon_plot,lat_plot,h_diff_CTRL');
    title('$\overline{\Delta h}$','interpreter','latex')
    format_fig(h,ax)
    ax = subplot(1,4,3);
    [~,h] = contourf(lon_plot,lat_plot,h_diff_prime');
    title('$\Delta h''$','interpreter','latex')
    format_fig(h,ax)
        ax = subplot(1,4,4);
    [~,h] = contourf(lon_plot,lat_plot,h_diff'-(h_diff_CTRL'+h_diff_prime'));
    title('$diff$','interpreter','latex')
    format_fig(h,ax)






function [] = format_fig(h,plt_num,max_val,min_val)

set(h,'edgecolor','none')
set(gca,'ydir','normal','fontsize',15)
colorbar
xlabel('deg')
ylabel('deg')

if nargin>2 % red white and blue colormap
    colormap(plt_num,rwb_map([max_val 0 min_val],100))
else
    colormap(plt_num,parula)
end


end