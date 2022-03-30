


str = 'v10_NCEP';
v_str = 'v';

years = 2004:2018;

for y = years
    nc_file = sprintf('~/Downloads/%s_%d.nc',str,y);
    ncinfo(nc_file)
    v = ncread(nc_file,v_str);
    v(1)
end










