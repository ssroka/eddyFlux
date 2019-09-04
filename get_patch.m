function [var_patch] = get_patch(var_str,srcD,srcJFM,norm_fac,patch_lon,patch_lat,lsm_patch,print_time)
    tic
    var_tot = cat(3,ncread(srcD,var_str),ncread(srcJFM,var_str)).*norm_fac; 
    var_patch = var_tot(patch_lon,patch_lat,:);
    if nargin>6
    var_patch(repmat(lsm_patch>0,1,1,size(var_patch,3)))  =NaN;
    end
    etime = toc;
    if print_time
        fprintf('Loaded %s in %2.1f sec\n',var_str,etime)
    end
end