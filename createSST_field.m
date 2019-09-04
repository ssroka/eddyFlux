function [sstf] = createSST_field(N,Lx,Ly,h)
%% User Input

x_range = [-10 10];
y_range = [-10 10];

% h = 0.1;

% characteristic eddy length scale
% Lx       = 4;
% Ly       = 4;

sst_hi = 300;
sst_lo = 295;

%% Begin

% ============== set SST ==========================================

% N = floor(diff(x_range)/h)+1;
dx = round(Lx/h);
dy = round(Ly/h);
% create eddy field
sst = zeros(N);

mu_x = 0;
mu_y = 0;
sig_x = 0.5;
sig_y = 0.7;

x = repmat(linspace(-Lx/2,Lx/2,dx),dy,1);
y = repmat([linspace(-Ly/2,Ly/2,dy)]',1,dx);
alpha = 0.1;
filter_flag = true;

[sst_mag,sstf_mag]=create_field(N,dx,dy,sig_x,sig_y,mu_x,mu_y,x,y,filter_flag,alpha);
sst  = (sst_mag-min(sst_mag(:)))*(sst_hi-sst_lo)+sst_lo; % 
sstf = (sstf_mag-min(sstf_mag(:)))*(sst_hi-sst_lo)+sst_lo; % 

% subplot(3,1,1)
% imagesc(sst)
% colorbar
% subplot(3,1,2)
% imagesc(sstf)
% colorbar
% subplot(3,1,3)
% imagesc(sst-sstf)
% colorbar


function [e,ef]=create_field(N,dx,dy,sig_x,sig_y,mu_x,mu_y,x,y,filter_flag,alpha)

e = zeros(N);

for i = 1:floor(N/dx)
    X = 1+(i-1)*dx:dx+(i-1)*dx;
    for j = 1:floor(N/dy)
        Y = 1+(j-1)*dy:dy+(j-1)*dy;
        sx = sig_x+rand*0.5;
        sy = sig_y+rand*0.5;
        mx = mu_x+rand*0.5;
        my = mu_y+rand*0.5;
        e(Y,X) = (-1).^(mod(i,2))*(-1).^(mod(j,2))*exp(-(x - mx).^2/2/sx^2-(y - my).^2/2/sy^2);
    end
end

if filter_flag
    % low pass filter
    ef = zeros(size(e));
    ef(:,1) = e(:,1);
    %  alpha = 0.1;
    for i = 2:size(e,1)
        for j = 2:size(e,2)
            ef(i,j) = 0.5*(ef(i-1,j) + alpha*(e(i,j)-ef(i-1,j)) +...
                ef(i,j-1) + alpha*(e(i,j)-ef(i,j-1)) );
        end
    end
end

end


end









