
clear;
clc;
close ;

for j = 5
    switch j
        case 1
            SSH_src = 'ssh_grids_v1812_2007122812.nc';
            ssha_str = 'SLA';
            time = ncread(SSH_src,'Time')*24*60*60;
            time = datetime(time,'ConvertFrom','epochtime','epoch','1985-01-01');
            
            subplot(1,2,j)
            ssha = ncread(SSH_src,ssha_str)*100;
            contourf(ssha(:,:,1),-40:10:80)
            colorbar
        case 2
            SSH_src = 'CCAR_recon_sea_level_20000103_20090627_v1.nc';
            ssha_str = 'ssha';
            time = ncread(SSH_src,'time')*24*60*60;
            time = datetime(time,'ConvertFrom','epochtime','epoch','1900-01-01');
            for i = 1:length(time)
                if (time(i).Year == 2007) && (time(i).Month == 12)  && (time(i).Day == 27)
                    subplot(1,2,j)
                    ssha = ncread(SSH_src,ssha_str);
                    contourf(ssha(:,:,1)',-40:10:80)
                    colorbar
                end
            end
        case 3
            % SSH_src = 'osu_cioss_weekly_msla_geovel_2007_v1.nc';
            % ssha_str = 'sea_surface_height_above_sea_level';
        case 4
%             SSH_src = 'global-reanalysis-phy-001-030-daily_1589857735177.nc';
            SSH_src = 'ssh_data_global/global-reanalysis-phy-001-030-daily_1589860194526.nc';
            ssha_str = 'zos';
            time = ncread(SSH_src,'time')*60*60;
            time = datetime(time,'ConvertFrom','epochtime','epoch','1950-01-01');
            lat = ncread(SSH_src,'latitude');
            lon = ncread(SSH_src,'longitude');
            SSH = ncread(SSH_src,ssha_str);
            save('ssh_data_global/global-reanalysis-phy-001-030-daily_20021201_20070401','SSH','lat','lon','time')
            contourf(SSH(:,:,1)')
            colorbar
        case 5 
            SSH_src = 'ssh_data_global/global-reanalysis-phy-001-030-daily_1592498841102.nc';
            ssha_str = 'zos';
            time_total = ncread(SSH_src,'time')*60*60;
            time_total = datetime(time_total,'ConvertFrom','epochtime','epoch','1950-01-01');
            lat = ncread(SSH_src,'latitude');
            lon = ncread(SSH_src,'longitude');
            SSH_total = ncread(SSH_src,ssha_str);
            % save in several batches
            inds = (time_total.Year<=2012) & ((time_total.Month==12) | (time_total.Month==1) | (time_total.Month==2) | (time_total.Month==3))  ;
            SSH = SSH_total(:,:,inds);
            time = time_total(inds);
            save('global-reanalysis-phy-001-030-daily_20071201_20120331','SSH','lat','lon','time')
            inds = (time_total.Year>=2012) & ((time_total.Month==12) | (time_total.Month==1) | (time_total.Month==2) | (time_total.Month==3))  ;
            SSH = SSH_total(:,:,inds);
            time = time_total(inds);
            save('global-reanalysis-phy-001-030-daily_20131201_20180331','SSH','lat','lon','time')
            contourf(SSH(:,:,1)')
            colorbar
            %{
            DOWNLOAD	NAME	DESCRIPTION	STANDARD NAME	UNITS
	
bottomT
Sea floor potential temperature
sea_water_potential_temperature_at_sea_floor
degrees_C
	
mlotst
Density ocean mixed layer thickness
ocean_mixed_layer_thickness_defined_by_sigma_theta
m
	
siconc
Ice concentration
sea_ice_area_fraction
1
	
sithick
Sea ice thickness
sea_ice_thickness
m
	
so
Salinity
sea_water_salinity
1e-3
	
thetao
Temperature
sea_water_potential_temperature
degrees_C
	
uo
Eastward velocity
eastward_sea_water_velocity
m s-1
	
usi
Sea ice eastward velocity
eastward_sea_ice_velocity
m s-1
	
vo
Northward velocity
northward_sea_water_velocity
m s-1
	
vsi
Sea ice northward velocity
northward_sea_ice_velocity
m s-1
	
zos
Sea surface height
sea_surface_height_above_geoid
m
            
            
            %}
    end
    
end




