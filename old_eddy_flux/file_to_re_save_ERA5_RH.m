ccc

filter_flag   = 'box'; % 'box' or 'zonal'

patch_str = 'Kur'; % 'GS'   'Kur'

% Kurishio
lat_bnds = [25 45];
lon_bnds = [130 170];
lat_sm_bnds = [30 40];

for year = 2003:2007
    
    filename_data =sprintf('/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/ERA5_patch_data_%d',year);
    tic
    load(filename_data)
    load_time = toc;
    
    rho_a = 1.2;     % kg m^-3
    c_p_air = 1000;  % J / kg / K
    Lv = 2.26e6;     % J/kg
    
    
    %{
    clear RH_patch qo_patch qa_patch
    
    qo_patch = SAM_qsatWater(SST_patch, P0_patch) ;
    mask = isnan(SST_patch);
    qo_patch(mask) = NaN;

    RH_patch = zeros(size(t2m_patch));
    qa_patch = zeros(size(t2m_patch));
    for tt = 1:size(RH_patch,3)
        T = t2m_patch(:,:,tt);
        TD = d2m_patch(:,:,tt);
        TC = T-273.15;
        TDC = TD-273.15;
        es = SAM_psatWater(t2m_patch(:,:,tt)); % Pa
        %         http://glossary.ametsoc.org/wiki/Clausius-clapeyron_equation
        %         es = 0.6112*exp((17.67.*TC)./(243.5+TC)); %kPa
        % Bolton, D. 1980. The computation of equivalent potential temperature.
        % Mon. Wea. Rev.. 108. 1046?1053.
        e = 0.6112*exp((17.67.*TDC)./(243.5+TDC))*1000; % Pa
        RH_patch(:,:,tt) = 100*e./es;
        
        r = 0.622 * e ./ (P0_patch(:,:,tt)-e);
        qa_patch(:,:,tt) = r./(1+r);
    end
    %}
    
    tic
    save(sprintf('/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/ERA5_patch_data_%d',year),...
        'Lv','rho_a','c_p_air','lat','lon','lat_bnds','lon_bnds','patch_lat','patch_lon',... % lat and lon vars
        'DT_patch','P0_patch','SST_patch','lsm_patch','slhf_patch','sshf_patch',... % ERA5 fields
        'u10_patch','v10_patch','time','t2m_patch','d2m_patch',...
        'U_mag','RH_patch','qa_patch','qo_patch')           % calculated fields
    save_time = toc;
    
    fprintf('saved year %d, load time: %4.2f, save time: %4.2f \n',year,load_time,save_time)
    
    
    
end