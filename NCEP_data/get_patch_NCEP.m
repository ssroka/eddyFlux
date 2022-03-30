function [var_patch] = get_patch_NCEP(var_str,src,patch_lon,patch_lat,lsm_patch,print_time)
    tic
    var_tot = double(ncread(src,var_str));
    var_patch = var_tot(patch_lon,patch_lat,:);
    var_patch(repmat(lsm_patch>0,1,1,size(var_patch,3)))=NaN;
    etime = toc;
    if print_time
        fprintf('Loaded %s in %2.1f sec\n',var_str,etime)
    end
end