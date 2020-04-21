
ccc

x = linspace(-5,5,50);
y = linspace(-5,5,50);

[X,Y] = meshgrid(x,y);

z = zeros(50,50);

z(1:4,1:4) = 1;

z(1:2,1:2) = NaN;


m = size(z,1);
n = size(z,1);
n_row = 3;
n_col = 3;
NaN_inds = isnan(z);

M = boxcar_filter_mat(m,n,n_row,n_col,NaN_inds);

[fz,~] = boxcar_filter(z,M);

subplot(121)
imagesc(x,y,z)
xlim([-5 -4])
ylim([-5 -4])
colorbar
subplot(122)
imagesc(x,y,fz)
xlim([-5 -4])
ylim([-5 -4])
colorbar











figure

x = linspace(-5,5,50);
y = linspace(-5,5,50);


[X,Y] = meshgrid(x,y);
sig = 1;

G = exp(-(X.^2+Y.^2)/(2*sig));

surf(x,y,G)
colorbar
colormap('bone')
set(gca,'GridLineStyle','none')
alpha 0.25

m = size(G,1);
n = size(G,1);
n_row = 10;
n_col = 10;
NaN_inds = isnan(G);

M = boxcar_filter_mat(m,n,n_row,n_col,NaN_inds);
fG = G;
for i = 1:5
[fG,~] = boxcar_filter(fG,M);
end
hold on

surf(x,y,reshape(fG,m,n))
colorbar
colormap('bone')
alpha 0.5
set(gca,'GridLineStyle','none')










