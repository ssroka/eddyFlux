clear;clc;close all

x = linspace(-5,5,10);

[X,Y] = meshgrid(x,x);

U = Y;
V = -X;



[cz,cav] = curl(X,Y,U,V);

subplot(1,2,1);

imagesc(x,x,cz); colorbar
hold on
q = quiver(X,Y,U,V);
title('\nabla x \tau')
set(gca,'fontsize',18)

[d] = divergence(X,Y,U,V);

subplot(1,2,2);
imagesc(x,x,d); colorbar
hold on
q = quiver(X,Y,U,V);
title('\nabla \cdot \tau')
set(gca,'fontsize',18)




