function [var_patch] = get_patch(var_str,srcD,DEC_ids,srcJFM,JFM_ids,norm_fac,patch_lon,patch_lat,lsm_patch,print_time)
    tic
    var_DEC = ncread(srcD,var_str);
    var_JFM = ncread(srcJFM,var_str);
    var_tot = cat(3,var_DEC(:,:,DEC_ids),var_JFM(:,:,JFM_ids)).*norm_fac; 
    var_patch = var_tot(patch_lon,patch_lat,:);
    if nargin>6
    var_patch(repmat(lsm_patch>0,1,1,size(var_patch,3)))  =NaN;
    end
    etime = toc;
    if print_time
        fprintf('Loaded %s in %2.1f sec\n',var_str,etime)
    end
end