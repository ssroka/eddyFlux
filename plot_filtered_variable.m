figure
subplot(4,2,[1 2])
contourf(lon_er,lat_er,A')
title('T_o')
set(gca,'ydir','normal')
set(gca,'fontsize',20)
colorbar

subplot(4,2,3)
contourf(lon_er,lat_er,(A_sm_x)')
title('x-dir filter T_o')

subplot(4,2,4)
contourf(lon_er,lat_er,(A-A_sm_x)')
title('T_o - x-dir filter T_o')

subplot(4,2,5)
contourf(lon_er,lat_er,(A_sm_y)')
title('y-dir filter T_o')

subplot(4,2,6)
contourf(lon_er,lat_er,(A-A_sm_y)')
title('T_o - y-dir filter T_o')

subplot(4,2,7)
contourf(lon_er,lat_er,(A_sm)')
title('0.5(y-dir T_o + x-dir T_o)')

subplot(4,2,8)
contourf(lon_er,lat_er,(A_prime)')
title('T_o - 0.5(y-dir T_o + x-dir T_o)')

for i = [3:8]
    subplot(4,2,i)
    set(gca,'ydir','normal')
    set(gca,'fontsize',20)
    colorbar
end

set(gcf,'position',[1     1   720   804],'color','w')

update_figure_paper_size()
print(sprintf('imgs/cmp_lanczos_%s_L_%d_%d',...
    'To',L/1000,year),'-dpdf')