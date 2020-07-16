ccc

% if you change the period length to anything that isn't a multiply of pi
% you get some leakage or ringing...the signal does not look as nice


cT = 700000; % cutoff period;

year = [2003];

box_limits = [34 41; 143 168];
box_limits = [30 42; 148 168];

data_src = '/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/';

%%

load(sprintf('%sERA5_patch_data_%d.mat',data_src,year));

lat = double(lat);
lon = double(lon);

lat_er = lat(patch_lat);
lon_er = lon(patch_lon);

lat_inds = lat_er>=box_limits(1,1) & lat_er<=box_limits(1,2);
lon_inds = lon_er>=box_limits(2,1) & lon_er<=box_limits(2,2);

lat_er = lat_er(lat_inds);
lon_er = lon_er(lon_inds);

m = size(SST_patch,1);
n = size(SST_patch,2);

d_lat = abs(lat(2)-lat(1));
d_lon = abs(lon(2)-lon(1));

m_per_deg = 111320; % only used to get the filter width;

dx = abs(d_lat*m_per_deg);
dy = abs(d_lon*m_per_deg);

cf = 1/cT;

s = nanmean(SST_patch(lon_inds,lat_inds,:),3)';

s = s - mean(s(:));

s = detrend(s); % detrend columns
s = detrend(s')'; % detrend rows

m = size(s,1);
n = size(s,2);

subplot(4,3,1)
contourf(lon_er,lat_er,s)
title(sprintf('detrend(ERA SST - mean(SST))'))

Fs = 1/dx; % samples per period, L samples for NT periods
Fy = [-m/2:m/2-1]/m*Fs;
Fx = [-n/2:n/2-1]/n*Fs;

[Xhat,Yhat] = meshgrid(Fx,Fy);

Zhat = fft2(s);

Zhat_shift = fftshift(Zhat);

Zhat_LP = Zhat_shift;
Zhat_HP = Zhat_shift;

inds_LP= sqrt(Xhat.^2+Yhat.^2)>cf;
Zhat_LP(inds_LP) = 0.0;

inds_HP = sqrt(Xhat.^2+Yhat.^2)<cf;
Zhat_HP(inds_HP) = 0.0;

Z_LP = ifft2(ifftshift(Zhat_LP));
Z_HP = ifft2(ifftshift(Zhat_HP));

Zhat_LP_shift = fftshift(fft2(Z_LP));
Zhat_HP_shift = fftshift(fft2(Z_HP));

Zhat_LP_HP = Zhat_LP_shift;
Zhat_LP_HP(inds_HP) = 0.0;

Zhat_HP_LP = Zhat_HP_shift;
Zhat_HP_LP(inds_LP) = 0.0;

Z_LP_HP = ifft2(ifftshift(Zhat_LP_HP));
Z_HP_LP = ifft2(ifftshift(Zhat_HP_LP));

s_prime = s-Z_LP;

s_p_hat = fftshift(fft2(s_prime));

s_p_hat_LP = s_p_hat;
s_p_hat_LP(inds_LP) = 0.0;

s_p_LP = ifft2(ifftshift(s_p_hat_LP));

Z_inv = ifft2(fft2(s));

subplot(4,3,2)
contourf(Fx,Fy,abs(Zhat_shift),'displayname','FFT(s)');
hold on 
title('FFT(SST)')

subplot(4,3,3)
contourf(Fx,Fy,abs(Zhat_LP),'--','displayname','FFT(s)');
title('low pass')

subplot(4,3,4)
contourf(Fx,Fy,abs(Zhat_HP),'--','displayname','FFT(s)');
title('high pass')

subplot(4,3,5)
contourf(lon_er,lat_er,abs(Z_LP))
title('low pass = bar(SST)')

subplot(4,3,6)
contourf(lon_er,lat_er,abs(Z_HP))
title('high pass')

subplot(4,3,7)
contourf(lon_er,lat_er,(Z_inv))
title('IFFT(FFT(SST))')

subplot(4,3,8)
contourf(lon_er,lat_er,abs(Z_LP_HP))
title('high pass(low pass)')

subplot(4,3,9)
contourf(lon_er,lat_er,abs(Z_HP_LP))
title('low pass(high pass)')

subplot(4,3,10)
contourf(lon_er,lat_er,s-abs(Z_LP))
title('SST'' = SST - bar(SST)')

subplot(4,3,11)
contourf(lon_er,lat_er,abs(s_p_LP))
title('low pass(SST'')')

subplot(4,3,12)
contourf(lon_er,lat_er,abs(Z_LP)+abs(Z_HP)-s)
title('low pass(SST'')')

for i = 1:12
    subplot(4,3,i)
    colorbar
end

set(gcf,'position',[440           1        1001         797])


update_figure_paper_size()

print(sprintf('/Users/ssroka/MIT/Research/eddyFlux/imgs/test_2D_fft'),'-dpdf')
