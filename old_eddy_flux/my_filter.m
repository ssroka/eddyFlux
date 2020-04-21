function [A_sm] = my_filter(A,M,Nx,Ny,N)

A_sm = A;
inds = isnan(A);
for j = 1:N
A_sm = ndnanfilter(A_sm,@rectwin,[Ny Nx]);
A_sm(inds) = NaN;
end



end





















