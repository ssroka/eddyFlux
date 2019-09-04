clear
close all
clc

% model = createpde();
% geometryFromEdges(model,@lshapeg);
% pdegplot(model,'EdgeLabels','on')
% ylim([-1.1,1.1])
% axis equal
% applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);
% specifyCoefficients(model,'m',0,...
%                           'd',0,...
%                           'c',1,...
%                           'a',0,...
%                           'f',1);
% generateMesh(model,'Hmax',0.25);
% results = solvepde(model);
% pdeplot(model,'XYData',results.NodalSolution)


model = createpde();
% Rectangle is code 3, 4 sides,
% followed by x-coordinates and then y-coordinates
geom = [3,4,0,1,1,0,0,0,1,1]';
g = decsg(geom);
geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
ylim([-1.1,1.1])
axis equal
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);
specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',@fcoeffunction);
generateMesh(model,'Hmax',0.25);
results = solvepde(model);
pdeplot(model,'XYData',results.NodalSolution)



function f = fcoeffunction(location,state)

% N = 3; % Number of equations
% nr = length(location.x); % Number of columns
% f = zeros(N,nr); % Allocate f
% 
% % Now the particular functional form of f
% f(1,:) = location.x - location.y + state.u(1,:);
% f(2,:) = 1 + tanh(state.ux(1,:)) + tanh(state.uy(3,:));
% f(3,:) = (5 + state.u(3,:)).*sqrt(location.x.^2 + location.y.^2);
Lx = 1; 
Ly = 1;
location
sst_hi = 300; % K
sst_lo = 295; % K
f = (sst_hi-sst_lo)*0.5*(cos(location.x*2*pi/Lx).*sin(location.y*2*pi/Ly)*2*pi/Lx)+...
    (sst_hi-sst_lo)*0.5*(sin(location.x*2*pi/Lx).*cos(location.y*2*pi/Ly)*2*pi/Ly);
% f = 1*ones(size(location.x));
end
