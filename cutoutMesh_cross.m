function [FEMesh, cutOutElements] = cutoutMesh_cross(L, Xk, Yk,maxMeshSize,minMeshSize)
gMatLength = 50; % Sufficiently Large

gd = zeros(gMatLength, 3);

%% Define the cross kirigami shape and the square substrate shape in 2D
c=1;
gd(1, c) = 2;
nedge=12;
gd(2, c) = nedge;
gd(3:2+nedge, c) = [-Xk,-Xk,Xk,Xk, Yk, Yk, Xk, Xk, -Xk, -Xk, -Yk,-Yk ];
gd(3+nedge:3+nedge*2-1, c) = [Xk, Yk, Yk, Xk,Xk,-Xk,-Xk,-Yk,-Yk,-Xk,  -Xk,Xk ];
dnn=20;
c=c+1;
gd(1:10, c) = [2,4, 0, (L/2)/dnn, 0, -(L/2)/dnn, (L/2)/dnn, 0, -(L/2)/dnn,0];
c=c+1;

HH=1.1 *(L/2);
gd(1:10, c) = [2,4, 0, HH, 0, -HH, HH, 0, -HH,0]; %% This is for rectangle

%% Define the model from the geometry
g = decsg(gd);
model = createpde;
geometryFromEdges(model,g);

%% Plot the geometry (for simple check)
figure(1)
pdegplot(model,'EdgeLabels','off')

%% Mesh the model

FEMesh = generateMesh(model,'Hmax', maxMeshSize, 'Hmin',minMeshSize);

kirigami_cutout=1:2; %Select the faces which form the kirigiami elements

figure(2)
pdeplot(model);
hold on
cutOutElements = findElements(FEMesh,'region','Face',kirigami_cutout); % Define the kirigiami elements
pdemesh(FEMesh.Nodes, FEMesh.Elements(:,cutOutElements),'EdgeColor','green');
hold off
axis equal
pp=1.1;
xlim([-pp*(L/2),pp*(L/2)])
ylim([-pp*(L/2),pp*(L/2)])
end