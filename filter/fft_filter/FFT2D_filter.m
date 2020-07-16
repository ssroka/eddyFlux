% =========================================================================
% 2D FFT Filter
% =========================================================================
% Created By    : Sydney Sroka
% Last Edited By: Sydney Sroka
% Date          : 07/14/2020
%
% ------------------------
% 
% This code applies the MATLAB built-in 2D FFT functions to a 2D field and
% returns the low pass and high pass fields per a user-supplied cutoff
% frequency.
% 
% INPUT:
%       s:          An mxn 2D array with uniform grid spacing (dx) in both
%                   directions.
% 
%       dx:         required when debug_flag is true. The spacing between the sample points.
% 
%       cf:         The cutoff frequency to apply (e.g. if you want to filter
%                   out features smaller than 100 dx in the low-pass field,
%                   then cf = 1/(100*dx) ).
% 
%       debug_flag: (optional) When true, creates contour plots of 
%                   intermediate output 
%  
%       x_vec:      (optional) This is a vector of n values evenly spaced by
%                    dx that is used in the debug plot
% 
%       y_vec:      (optional) This is a vector of m values evenly spaced by
%                    dx that is used in the debug plot
% 
% OUTPUT:
%      S_LP:        The low-pass filtered field.
% 
%      S_HP:        The high-pass filtered field.
% 
% 


function [Z_LP,Z_HP] = FFT2D_filter(s,dx,cf,debug_flag,x_vec,y_vec)

if any(isnan(s(:)))
    error('cannot process fft2 with NaN in input field')
end

m = size(s,1);
n = size(s,2);

Fs = 1/dx; % sampling frequency
Fy = [-m/2:m/2-1]/m*Fs+mod(m,2)*Fs/m/2;
Fx = [-n/2:n/2-1]/n*Fs+mod(n,2)*Fs/n/2;

[Xhat,Yhat] = meshgrid(Fx,Fy);

% take the N-point 2D FFT with the MATLAB built-in
Zhat = fft2(s);

% shift the output so that the zero freuency content is centered
Zhat_shift = fftshift(Zhat);

% ---------- low pass --------------------

% initialize the LP filtered field
Zhat_LP = Zhat_shift;
% indentify the indices outside of the circle centered at the origin with
% radius cf, that represent frequencies higher than the cutoff
inds_LP= sqrt(Xhat.^2+Yhat.^2)>cf;
% remove the high-frequency content
Zhat_LP(inds_LP) = 0.0;

% inverse shift and inverse fft to recover the low-pass field
Z_LP = ifft2(ifftshift(Zhat_LP));

% ---------- high pass --------------------

% initialize the HP filtered field
Zhat_HP = Zhat_shift;
% indentify the indices inside of the circle centered at the origin with
% radius cf, that represent frequencies lower than the cutoff
inds_HP = sqrt(Xhat.^2+Yhat.^2)<cf;
% remove the low-frequency content
Zhat_HP(inds_HP) = 0.0;

% inverse shift and inverse fft to recover the high-pass field
Z_HP = ifft2(ifftshift(Zhat_HP));


if (nargin > 3) &&  debug_flag
    
    if nargin<4
        x_vec = (0:n-1)*dx;
    end
    if nargin<5
        y_vec = (0:m-1)*dx;
    end
    
    % FFT and FFT shift the high pass and low pass fields
    Zhat_LP_shift = fftshift(fft2(Z_LP));
    Zhat_HP_shift = fftshift(fft2(Z_HP));
    
    % high pass filter the low pass field
    Zhat_LP_HP = Zhat_LP_shift;
    Zhat_LP_HP(inds_HP) = 0.0;
    
    % low pass filter the low pass field
    Zhat_LP_LP = Zhat_LP_shift;
    Zhat_LP_LP(inds_LP) = 0.0;
    
    % low pass filter the high pass field
    Zhat_HP_LP = Zhat_HP_shift;
    Zhat_HP_LP(inds_LP) = 0.0;
    
    Z_LP_HP = ifft2(ifftshift(Zhat_LP_HP));
    Z_LP_LP = ifft2(ifftshift(Zhat_LP_LP));
    Z_HP_LP = ifft2(ifftshift(Zhat_HP_LP));
    
    s_prime = s-Z_LP;
    
    s_p_hat = fftshift(fft2(s_prime));
    
    s_p_hat_LP = s_p_hat;
    s_p_hat_LP(inds_LP) = 0.0;
    
    s_p_LP = ifft2(ifftshift(s_p_hat_LP));
    
    Z_inv = ifft2(fft2(s));
    
    subplot(4,3,1)
    contourf(x_vec,y_vec,s)
    title('input field','interpreter','latex')
    
    subplot(4,3,2)
    contourf(Fx,Fy,abs(Zhat_shift),'displayname','FFT(s)');
    hold on
    % plot circle at cutoff
    th = linspace(0,2*pi);
    plot(cf*cos(th),cf*sin(th),'r','linewidth',2);
    xlabel(' -1/(2*dx) to 1/(2*dx) ')
    ylabel(' -1/(2*dx) to 1/(2*dx) ')
    hold on
    title('FFT, red line shows $cf$','interpreter','latex')
    
    subplot(4,3,3)
    contourf(Fx,Fy,abs(Zhat_LP));
        xlabel(' -1/(2*dx) to 1/(2*dx) ')
    ylabel(' -1/(2*dx) to 1/(2*dx) ')
    title('LP: FFT(s) with high frequencies removed ','interpreter','latex')
    
    subplot(4,3,4)
    contourf(Fx,Fy,abs(Zhat_HP));
        xlabel(' -1/(2*dx) to 1/(2*dx) ')
    ylabel(' -1/(2*dx) to 1/(2*dx) ')
    title('HP: FFT(s) with low frequencies removed ','interpreter','latex')
    
    subplot(4,3,5)
    contourf(x_vec,y_vec,abs(Z_LP))
    title('low pass output = $\mathcal{F}_{LP}(s)$','interpreter','latex')
    
    subplot(4,3,6)
    contourf(x_vec,y_vec,abs(Z_HP))
    title('high pass output = $\mathcal{F}_{HP}(s)$','interpreter','latex')
    
    subplot(4,3,7)
    contourf(x_vec,y_vec,abs(Z_LP_HP))
    title('$\mathcal{F}_{HP}(\mathcal{F}_{LP}(s))$ (should vanish)','interpreter','latex')
    
    subplot(4,3,8)
    contourf(x_vec,y_vec,abs(Z_HP_LP))
    title('$\mathcal{F}_{LP}(\mathcal{F}_{HP}(s))$ (should vanish)','interpreter','latex')
    
    subplot(4,3,9)
    contourf(x_vec,y_vec,s-abs(Z_LP))
    title('perturbation to the mean s'' = s - $\mathcal{F}_{LP}(s)$','interpreter','latex')
    
    subplot(4,3,10)
    contourf(x_vec,y_vec,abs(s_p_LP))
    title('$\mathcal{F}_{LP}(s'')$ (should vanish)','interpreter','latex')
    
    subplot(4,3,11)
    contourf(x_vec,y_vec,abs(Z_LP_LP))
    title('$\mathcal{F}_{LP}(\mathcal{F}_{LP}(s))$','interpreter','latex')
    
    subplot(4,3,12)
    contourf(x_vec,y_vec,abs(Z_LP_LP) - abs(Z_LP))
    title('$\mathcal{F}_{LP}(\mathcal{F}_{LP}(s)) - \mathcal{F}_{LP}(s)$ (should vanish)','interpreter','latex')
    
    for i = 1:12
        subplot(4,3,i)
        colorbar
    end
end



end
