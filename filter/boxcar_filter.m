function  [A_sm,A_prime] = boxcar_filter(A,M)
INDS = isnan(A);
num_nan = sum(INDS,2);
orig_mean_num = diag(M);
A(INDS) = 0;
A_sm = reshape(M*A(:),size(A,1),size(A,2));
A_prime = reshape(A(:) - M*A(:),size(A,1),size(A,2));
A_sm(INDS) = NaN;
A_prime(INDS) = NaN;



