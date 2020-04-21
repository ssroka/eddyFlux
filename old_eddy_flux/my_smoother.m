function [A_mean, A_prime] = my_smoother(A,n_row,n_col)

% make sure row/col length is odd
if mod(n_col,2)==0
    n_col = n_col + 1;
end
if mod(n_row,2)==0
    n_row = n_row + 1;
end

LR = (n_col-1)/2;
UD = (n_row-1)/2;

m = size(A,1);
n = size(A,2);

A_prime = zeros(size(A));
A_mean= zeros(size(A));

for i = 1:m
    st_row = max(1,i-UD);
    end_row = min(m,i+UD);
    for j = 1:n
        st_col = max(1,j-LR);
        end_col = min(n,j+LR);
        A_sub = A(st_row:end_row,st_col:end_col);
        A_mean(i,j) = nanmean(A_sub(:));
        A_prime(i,j) = A(i,j) - A_mean(i,j);
    end
end


