clear;close all;clc

test_mat = ones(10,11);
test_mat(5,5) = 2;
test_mat(1,1) = NaN;
test_mat(2,1) = NaN;
test_mat(6,6) = NaN;

imagesc(test_mat)

Nx = 2;
Ny = 2;

m = size(test_mat,1);
n = size(test_mat,2);

NaN_inds = isnan(test_mat(:));

[M] = my_smoother_mat(m,n,Ny,Nx,NaN_inds);
[sm_mat1,p_mat1] = smooth_mat(test_mat,M);
[sm_mat2] = ndnanfilter(test_mat,[],[2 2]);
p_mat2 = test_mat - sm_mat2;

[sm_pmat1,pp_mat1] = smooth_mat(p_mat1,M);
[sm_pmat2] = ndnanfilter(p_mat2,[],[2 2]);
pp_mat2 = p_mat2 - sm_pmat2;

[sm_sm_mat1,~] = smooth_mat(sm_mat1,M);
[sm_sm_mat2] = ndnanfilter(sm_mat2,[],[2 2]);


subplot(3,4,1)
imagesc(test_mat)
colorbar
title('raw')
subplot(3,4,2)
imagesc(sm_mat1)
colorbar
title('smoothed with M')
subplot(3,4,3)
imagesc(sm_mat2)
colorbar
title('smoothed without M')
subplot(3,4,4)
imagesc(sm_mat2-sm_mat1)
colorbar
title('smoothed diff')
subplot(3,4,5)
imagesc(test_mat)
colorbar
title('raw')
subplot(3,4,6)
imagesc(p_mat1)
colorbar
title('prime with M')
subplot(3,4,7)
imagesc(p_mat2)
colorbar
title('prime without M')
subplot(3,4,8)
imagesc(p_mat2-p_mat1)
colorbar
title('prime diff')
subplot(3,4,9)
imagesc(sm_pmat1)
colorbar
title('sm pmat 1')
subplot(3,4,10)
imagesc(sm_pmat2)
colorbar
title('sm pmat 2')
subplot(3,4,11)
imagesc(sm_sm_mat1-sm_mat1)
colorbar
title('diff sm sm mat 1 with sm mat 1')
subplot(3,4,12)
imagesc(sm_sm_mat2-sm_mat2)
colorbar
title('diff sm sm mat 2 with sm mat 2')







