function [A_sm,A_prime,A_sm_x,A_sm_y] = lanczos_filter(A,dx,cf)

A_sm_x = zeros(size(A,1),size(A,2));
for i = 1:size(A,1)
    nan_inds = isnan(A(i,:));
    A_sm_x(i,:) = lanczosfilter(A(i,:),dx,cf);
    A_sm_x(i,nan_inds) = NaN;
end

% A_sm = A_sm_x;

A_sm_y = zeros(size(A,1),size(A,2));
for i = 1:size(A,2)
    nan_inds = isnan(A(:,i)');
    A_sm_y(:,i) = lanczosfilter(A(:,i)',dx,cf)';
    A_sm_y(nan_inds,i) = NaN;
end

A_sm = 0.5*(A_sm_x+A_sm_y);

A_prime = A - A_sm;

end

































