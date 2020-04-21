clear
close all
clc


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
                          'f',@f_pot);
generateMesh(model,'Hmax',0.05);
results = solvepde(model);
% pdeplot(model,'XYData',results.NodalSolution)
pdeplot(model,'XYData',results.XGradients)
% plot3(results.Mesh.Nodes(1,:),results.Mesh.Nodes(2,:),results.NodalSolution,'o')
xn = results.Mesh.Nodes(1,:)';
yn = results.Mesh.Nodes(2,:)';
u = xn.^2.*yn.*(xn-1).*(yn-1);
v = yn.^2.*xn.*(xn-1).*(yn-1);
figure
plot3(xn,yn,u,'bo')
hold on
plot3(xn,yn,v,'ro')






return
model2 = createpde();
% Rectangle is code 3, 4 sides,
% followed by x-coordinates and then y-coordinates
geom = [3,4,0,1,1,0,0,0,1,1]';
g = decsg(geom);
geometryFromEdges(model2,g);
pdegplot(model2,'EdgeLabels','on')
ylim([-1.1,1.1])
axis equal
applyBoundaryCondition(model2,'dirichlet','Edge',1:model2.Geometry.NumEdges,'u',0);
specifyCoefficients(model2,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',@f_sf);
generateMesh(model2,'Hmax',0.05);
results2 = solvepde(model2);
pdeplot(model2,'XYData',results2.YGradients+results.XGradients)
figure
xn = results.Mesh.Nodes(1,:);
yn = results.Mesh.Nodes(2,:);
u = xn.^2.*yn.*(xn-1).*(yn-1);
x = linspace(0,1);
[X,Y] = meshgrid(x,x);
F = scatteredInterpolant(xn',yn',results.NodalSolution);
Fv = F(X,Y);
[Fx,Fy] = gradient(Fv,x(2)-x(1),x(2)-x(1));
imagesc(Fx)
plot3(xn,yn,u,'ro')
figure
xn = results.Mesh.Nodes(1,:);
yn = results.Mesh.Nodes(2,:);
plot3(xn,yn,results2.YGradients-results.XGradients,'o')
plot3(results.Mesh.Nodes(1,:),results.Mesh.Nodes(2,:),results.NodalSolution,'o')



function f = f_pot(location,state)
x = location.x;
y = location.y;
f = -(6.*x.^2.*y.^2+4.*x.*y-5.*x.^2.*y-5.*x.*y.^2);
end

function f = f_sf(location,state)
x = location.x;
y = location.y;
f = 2.*x.*y.^3-y.^3-2.*y.^2.*x+y.^2-(2.*x.^3.*y-x.^3-2.*x.^2.*y+x.^2);
end