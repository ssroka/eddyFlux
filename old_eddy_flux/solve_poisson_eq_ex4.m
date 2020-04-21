clear
close all
clc

%{
V = < y , 3 x >
%}


model = createpde();
% Rectangle is code 3, 4 sides,
% followed by x-coordinates and then y-coordinates
geom = [3,4,0,1,1,0,0,0,1,1]';
g = decsg(geom);
geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
ylim([-0.1,1.1])
axis equal

% bc eqs
e1_bc = @(location,state) 0.5.*(location.y.^2-3.*location.x.^2);
e2_bc = @(location,state) 0.5.*(1+location.y.^2);
e3_bc = @(location,state) 0.5.*(location.x.^2+1);
e4_bc = @(location,state) 0.5.*(location.y.^2);



applyBoundaryCondition(model,'dirichlet','Edge',1,'u',e1_bc);
applyBoundaryCondition(model,'dirichlet','Edge',2,'u',e1_bc);
applyBoundaryCondition(model,'dirichlet','Edge',3,'u',e1_bc);
applyBoundaryCondition(model,'dirichlet','Edge',4,'u',e1_bc);
specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',@f_pot);
generateMesh(model,'Hmax',0.05);
results = solvepde(model);
figure
% pdeplot(model,'XYData',results.NodalSolution)
subplot(2,1,1)
pdeplot(model,'XYData',results.YGradients)
title('u')
subplot(2,1,2)
pdeplot(model,'XYData',-results.XGradients)
title('v')


function f = f_pot(location,state)
f = 2*ones(size(location.x));
end

