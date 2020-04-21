clear;close all;clc
addpath('~/Documents/MATLAB/util/ndnanfilter/')
% I = imread('cameraman.tif');
load('2003_model_eddy_min_no_eddy_results','DT_eddy','patch_lat','patch_lon','lat','lon')
filt_L = 3;

Nx = 19;
Ny = 19;

[M] = my_smoother_mat(sum(patch_lon),sum(patch_lat),Ny,Nx);

I = DT_eddy;
inds = isnan(I);
filt_runs = [1 5 10 25 50];

% [Y] = ndnanfilter(I,[],[filt_L filt_L]);
Y = smooth_mat(I,M);

% I_prime_blur1 = ndnanfilter(I-Y,[],[filt_L filt_L]);
I_prime_blur1 = smooth_mat(I-Y,M);

I_prime_blur = I_prime_blur1;
I_prime_blur1(inds) = NaN;

for j = 1:length(filt_runs)
    
% Iblur = imgaussfilt(I,filt_L);
% Iblur_1 = imgaussfilt(Iblur,filt_L);

% Iblur = ndnanfilter(I,@rectwin,[filt_L filt_L]);
Iblur = smooth_mat(I,M);
Iblur(inds) = NaN;
% Iblur_1 = ndnanfilter(Iblur,[],[filt_L filt_L]);
Iblur_1 = smooth_mat(Iblur,M);
Iblur_1(inds) = NaN;

for i = 1:filt_runs(j)
% Iblur = imgaussfilt(Iblur,2);
% Iblur_1 = imgaussfilt(Iblur_1,2);

% I_prime_blur = ndnanfilter(I_prime_blur,[],[filt_L filt_L]);
I_prime_blur = smooth_mat(I_prime_blur,M);
I_prime_blur(inds) = NaN;

% Iblur = ndnanfilter(Iblur,[],[filt_L filt_L]);
Iblur = smooth_mat(Iblur,M);
Iblur(inds) = NaN;
% Iblur_1 = ndnanfilter(Iblur_1,[],[filt_L filt_L]);
Iblur_1 = smooth_mat(Iblur_1,M);
Iblur_1(inds) = NaN;
end
max_er(j) = max(abs(Iblur(:)-Iblur_1(:)));
end

figure

plot(filt_runs,max_er,'-o')
xlabel('filter runs')
ylabel('max diff btwn successive fields')

figure

subplot(3,3,1)
contourf(lon(patch_lon),lat(patch_lat),I')
title('$$\phi$$','interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,2)
contourf(lon(patch_lon),lat(patch_lat),Y')
title('$$<\phi>$$','interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,3)
contourf(lon(patch_lon),lat(patch_lat),I'-Y')
title('$$\phi'' = \phi - <\phi>$$','interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,4)
contourf(lon(patch_lon),lat(patch_lat),Iblur_1')
title(['$$<...<\phi>...>$$' sprintf(', %d filters',filt_runs(end))],'interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,5)
contourf(lon(patch_lon),lat(patch_lat),I'-Iblur_1')
title('$$\phi - <...<\phi>...>$$','interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,6)
contourf(lon(patch_lon),lat(patch_lat),I_prime_blur1')
title(['$$<\phi''>$$' sprintf(', %d filters',filt_runs(end))],'interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,7)
contourf(lon(patch_lon),lat(patch_lat),I_prime_blur')
title(['$$<...<\phi''>...>$$' sprintf(', %d filters',filt_runs(end))],'interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,8)
contourf(lon(patch_lon),lat(patch_lat),Y' - Iblur_1')
title(['$$<\phi> - <...<\phi>...>$$' sprintf(', %d filters',filt_runs(end))],'interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')
subplot(3,3,9)
contourf(lon(patch_lon),lat(patch_lat),Iblur_1'-Iblur')
title([ sprintf(' %d filters',filt_runs(end)+1) ' - '...
        sprintf(' %d filters',filt_runs(end))],'interpreter','latex')
colorbar;set(gca,'fontsize',20)
xlabel('lon deg')
ylabel('lat deg')


set(gcf,'position',[44         112        1124         693],'color','w')
set(gcf,'numbertitle','off','name','2003 DT')


