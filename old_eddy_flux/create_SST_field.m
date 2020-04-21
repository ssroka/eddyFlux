function [SST_field,SST_CTRL] = create_SST_field(ctrl_flag, SST_bar,D_SST,y_sst_km,l_sst_km,l_an_km,L_km,eddy_size,X,Y,n_eddies,rand_seed,filter_flag,L)
SST_CTRL_0 = SST_bar-D_SST/2*tanh((Y-y_sst_km)/l_sst_km);
if exist('rand_seed','var')
    rng(rand_seed)
end
if ~ctrl_flag
    SST_eddies_T = zeros(size(X));
    for i = 1:n_eddies
        x_c= rand*L_km;
        y_c= randn*y_sst_km/3+y_sst_km;
        sig_c = rand*diff(eddy_size)+eddy_size(1);
        SST_eddies_T = SST_eddies_T + exp((-(X-x_c).^2-(Y-y_c).^2)/2/sig_c^2);
    end
    
    % remove zonal mean
    SST_eddies_T = SST_eddies_T - repmat(mean(SST_eddies_T,2),1,size(X,2));
    a = fzero(@(a) std(reshape(SST_eddies_T,numel(SST_eddies_T),1)*a)-3,4);
    SST_eddies_T = SST_eddies_T*a;
    G = exp(-(Y-y_sst_km).^2./l_an_km.^2);
    SST_field = SST_CTRL_0 + SST_eddies_T.*G;
else
    SST_field = SST_CTRL_0;
end

if strcmp('box',filter_flag)
    X_dist = abs(X(end,end))-(X(end,1));
    Y_dist = abs(Y(1,1)-Y(end,1));
    
    Nx     = floor(X_dist/L);
    Ny     = floor(Y_dist/L);
    
    SST_CTRL = smooth2a(SST_field,Ny, Nx);
else
    SST_CTRL = SST_CTRL_0;
end





end