function [U,V,P] = create_baroclinic_wave(P0,delta_P,centerline_phi,L_km,N,lambda,rho_a)

press_AVG = P0;
press_HI = + delta_P; 
press_LO = - delta_P;  

% create two Gaussians 
sig = lambda/3;

x = linspace(-L_km/2,L_km/2,N);
[X, Y] = meshgrid(x,x);


loc_1 = [-lambda/2 0];
loc_2 = [lambda/2 0];

P_mean = press_AVG*ones(size(X));

t = press_LO*exp( (-(X-loc_1(1)).^2-(Y-loc_1(2)).^2)/(2*sig^2) );
t = t + press_HI*exp( (-(X-loc_2(1)).^2-(Y-loc_2(2)).^2)/(2*sig^2) );
P = P_mean + t;


Omega = 7.2921e-5; % rad/s
R = 6373; % radius of the earth in km
d_theta = L_km/R;

[dPdx,dPdy] = gradient(P,(x(2)-x(1))*1000); % gradient in Pa/m

phi = Y/R+(centerline_phi/180*pi);

f = 2*Omega.*sin(phi);


V = 1./rho_a./f.*dPdx;
U = -1./rho_a./f.*dPdy;

end







