function [M] = my_smoother_mat(m,n,n_row,n_col,NaN_inds)

% make sure row/col length is odd
if mod(n_col,2)==0
    n_col = n_col + 1;
end
if mod(n_row,2)==0
    n_row = n_row + 1;
end

LR = (n_col-1)/2;
UD = (n_row-1)/2;

IND = false(m*n,m*n);
for j = 1:n
    st_col = max(1,j-LR);
    end_col = min(n,j+LR);
    for i = 1:m
        st_row = max(1,i-UD);
        end_row = min(m,i+UD);
        A = false(m,n);
        A(st_row:end_row,st_col:end_col) =true;
        A(NaN_inds) = false;
        IND(i+(j-1)*m,:) = A(:)'&~NaN_inds(:)';
       
    end
end

M = 1./repmat(sum(IND,2),1,n*m).*IND;




