% ccc

figure

rng(0)
h_bar = 10*rand(5,5,5);
h_prime = rand(5,5,5);

U_bar = 10*rand(5,5,5);
To_prime= 2*repmat(rand(5,5),1,1,5);

b = 1.5;

q1 = U_bar.*h_bar;
q2 = U_bar.*h_prime;
q3 = b*To_prime.*h_bar;
q4 = b*To_prime.*h_prime;


r1 = (mean(q3,3)+mean(q4,3))./(mean(q1,3)+mean(q2,3));
r2 = mean(b*To_prime./U_bar,3);
contour(r1./r2)
colorbar
set(gcf,'Position',[ 48   332   560   420])














