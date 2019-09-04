function [f_bar,f_prime] = calc_anomaly(f,L)
nr = size(f,1);
nc = size(f,2);

field_borders = [f(end-L:end,end-L:end) f(end-L:end,:) f(end-L:end,1:L+1);
                 f(:,end-L:end) f f(:,1:L+1);
                 f(1:L+1,end-L:end) f(1:L+1,:) f(1:L+1,1:L+1)];

[X, Y] = meshgrid(1:nc+2*(L+1),1:nr+2*(L+1));
f_bar   = zeros(nr,nc);
f_prime = zeros(nr,nc);
% set up matrix multiplication
tic
parfor ii = 1:nr*nc
    [row_i,col_i] = ind2sub([nr nc],ii);
    row_i = row_i + L+1;
    col_i = col_i + L+1;
    get_neighbour_ids = (X-X(row_i,col_i)).^2 + (Y-Y(row_i,col_i)).^2 < L.^2;
    f_bar(ii)   = mean(field_borders(get_neighbour_ids));
    f_prime(ii) = field_borders(row_i,col_i)-mean(field_borders(get_neighbour_ids));
end
toc
disp('done')








end