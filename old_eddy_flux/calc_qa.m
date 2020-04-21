function [qa] = calc_qa(sst,p_surf,p_aloft,DT,RH)
% p_surf in Pa
% DT in K
% sst in K
% RH [%]

T_surf = sst-DT;

e_sat = SAM_psatWater(T_surf);
e = RH/100*e_sat;
r = 0.622 * e ./ max(e, p_surf-e);
qa = r./(1+r);

end

% rho_dry = 1.22; %kg m^-3
% R_w =  461.52; % J /kg /K
% R_dry = 287;   % J /kg /K

% e = rho_v*R_w*T_surf;
% p_dry_air = rho_dry*R_dry*T_surf;

% L_v = 2260000; % J/ kg
% T_aloft = 1./(R_w/L_v*log(p_surf/p_aloft)+1./T_surf);

% page 78 of P&K
