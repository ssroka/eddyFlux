
clear;
close all
clc


syms a CD rho_a real
syms T u v real
syms dTdx dTdy dudx dudy dvdx dvdy real

UX = 0.5*(u^4+u^2*v^2)^(-0.5)*(4*u^3*dudx+2*u*dudx*v^2+2*v*dvdx*u^2);
UY = 0.5*(u^4+u^2*v^2)^(-0.5)*(4*u^3*dudy+2*u*dudy*v^2+2*v*dvdy*u^2);

VX = 0.5*(v^4+v^2*u^2)^(-0.5)*(4*v^3*dvdx+2*v*dvdx*u^2+2*u*dudx*v^2);
VY = 0.5*(v^4+v^2*u^2)^(-0.5)*(4*v^3*dvdy+2*v*dvdy*u^2+2*u*dudy*v^2);


D = rho_a*CD*(a*dTdx*sqrt(u^4+u^2*v^2)+(1+a*T)*UX+...
              a*dTdy*sqrt(v^4+v^2*u^2)+(1+a*T)*VY);

C = rho_a*CD*(-a*dTdx*sqrt(v^4+v^2*u^2)-(1+a*T)*VX+...
              a*dTdy*sqrt(u^4+u^2*v^2)+(1+a*T)*UY);
          
pretty(D)

AW = (dTdx*u/(sqrt(u^2+v^2))+dTdy*v/(sqrt(u^2+v^2)))*...
      [u/(sqrt(u^2+v^2)),v/(sqrt(u^2+v^2))];
  
CW = [dTdx,dTdy] - AW;

pretty(AW)


pretty(CW)