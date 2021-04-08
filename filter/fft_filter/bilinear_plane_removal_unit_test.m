
ccc;

Nx = 10;
Ny = 20;

[X, Y] = meshgrid(-Nx/2:Nx/2,-Ny/2:Ny/2);

n = size(X,1);
m = size(X,2);

s = 4*X.^2+6*X+5*Y+3;

subplot(1,2,1)
surf(X,Y,s)
hold on
% use Least Squares to find the bilinear plane

DesMat = [X(:) Y(:) ones(numel(X),1)];

betas = DesMat\s(:)

s_LS = reshape(DesMat*betas,n,m);
surf(X,Y,s_LS)

subplot(1,2,2)
surf(X,Y,s-s_LS)
















