function [N,M] = get_num_subplots(P)
% P is the number of plots to make




N = floor(sqrt(P));
M = floor(sqrt(P));

row_inc = true;

while N*M < P
    
    if row_inc
        M = M + 1;
        row_inc = false;
    else
        N = N + 1;
        row_inc = true;
    end
    
end











end
