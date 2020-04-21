clear
close all
clc

% Rectangle is code 3, 4 sides,
% followed by x-coordinates and then y-coordinates
geom = [3,4,0,1,1,0,0,0,1,1]';
g = decsg(geom);

%% pot

model = createpde();

geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
ylim([-0.1,1.1])
axis equal

% bc eqs
e1_bc = @(location,state) location.y+2*location.x;
e2_bc = @(location,state) location.y+3*location.x;
applyBoundaryCondition(model,'neumann','Edge',1,'g',e1_bc);
applyBoundaryCondition(model,'neumann','Edge',3,'g',e1_bc);
applyBoundaryCondition(model,'neumann','Edge',2,'g',e2_bc);
applyBoundaryCondition(model,'neumann','Edge',4,'g',e2_bc);
specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',@f_pot);
generateMesh(model,'Hmax',0.05);
results = solvepde(model);
subplot(3,1,1)
% pdeplot(model,'XYData',results.NodalSolution)
pdeplot(model,'XYData',results.XGradients)

%% SF
return
model2 = createpde();

geometryFromEdges(model2,g);
% pdegplot(model2,'EdgeLabels','on')
% ylim([-0.1,1.1])
% axis equal

% bc eqs
e2_bc = @(location,state) 0.5.*(location.y.^2-3.*location.x.^2);
applyBoundaryCondition(model2,'dirichlet','Edge',1:4,'u',e2_bc);
specifyCoefficients(model2,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',@f_sf);
generateMesh(model2,'Hmax',0.05);
results2 = solvepde(model2);
subplot(3,1,2)
% pdeplot(model,'XYData',results.NodalSolution)
pdeplot(model2,'XYData',results2.XGradients)






function f = f_pot(location,state)
f = -3*ones(size(location.x));
end
function f = f_sf(location,state)
f = 2*ones(size(location.x));
end

