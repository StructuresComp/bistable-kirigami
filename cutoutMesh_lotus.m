function [FEMesh, cutOutElements] = cutoutMesh_lotus(innerRadius, outerRadius, ...
    Nstrips, dTheta, ...
    maxMeshSize, minMeshSize)

% innerRadius = inner radius of the cutout circle
% outerRadius = outer radius of the cutout circle
% Nstrips = number of strips/cuts
% dTheta = radial angle of each cut
% maxMeshSize = max size of an element
% minMeshSize = min size of an element

theta = linspace(0, 2*pi, Nstrips+1);
theta = theta(1:end-1); % 0 and 2*pi overlaps

%% Figure out the points on the quadrilateral cuts

% Additional points to capture the circular nature of the circle
Ndiscrete = 4;

ptQuads_x = zeros(Nstrips, (2*Ndiscrete));
ptQuads_y = zeros(Nstrips, (2*Ndiscrete));

for c=1:Nstrips
    th_1 = theta(c) - dTheta/2;
    th_2 = theta(c) + dTheta/2;
    
    th_arr = linspace(th_1, th_2, Ndiscrete);
    innerEdge_x = innerRadius * cos(th_arr);
    innerEdge_y = innerRadius * sin(th_arr);

    th_arr = linspace(th_2, th_1, Ndiscrete);
    outerEdge_x = outerRadius * cos(th_arr);
    outerEdge_y = outerRadius * sin(th_arr);
    
    % Store
    ptQuads_x(c, 1:Ndiscrete) = innerEdge_x;
    ptQuads_x(c, Ndiscrete+1:2*Ndiscrete) = outerEdge_x;
    ptQuads_y(c, 1:Ndiscrete) = innerEdge_y;
    ptQuads_y(c, Ndiscrete+1:2*Ndiscrete) = outerEdge_y;
end

%% Figure out the points on the circle minus quads
% Create a polygon
ptPoly_x = zeros(Nstrips*(2*Ndiscrete), 1);
ptPoly_y = zeros(Nstrips*(2*Ndiscrete), 1);

for c=1:Nstrips
    
    c2 = c + 1;
    if c==Nstrips
        c2 = 1;
    end
    
    th_1 = theta(c)  + dTheta/2;
    th_2 = theta(c2) - dTheta/2;
    if (th_2 < 0)
        th_2 = th_2 + 2*pi;
    end
    
    th_arr = linspace(th_1, th_2, Ndiscrete);
    outerEdge_x = outerRadius * cos(th_arr);
    outerEdge_y = outerRadius * sin(th_arr);
    
    % Store
    counter_i = (c-1)*(2*Ndiscrete) + 1;
    counter_f = (c-1)*(2*Ndiscrete) + Ndiscrete;
    ptPoly_x(counter_i:counter_f) = ptQuads_x(c, 1:Ndiscrete);
    ptPoly_y(counter_i:counter_f) = ptQuads_y(c, 1:Ndiscrete);

    ptPoly_x(counter_f+1:counter_f+Ndiscrete) = outerEdge_x;
    ptPoly_y(counter_f+1:counter_f+Ndiscrete) = outerEdge_y;
end

gMatLength = 1 + 1 + numel(ptPoly_x) * 2;
gd = zeros(gMatLength, Nstrips+1);
% Embed the quadrilateral cuts in the geometry
for c=1:Nstrips
    gd(1, c) = 2;
    gd(2, c) = numel(ptQuads_x(c,:));
    
    counter_i = 3;
    counter_f = 2 + numel(ptQuads_x(c,:));
    gd(counter_i:counter_f, c) = ptQuads_x(c,:);

    counter_i = counter_f + 1;
    counter_f = (counter_i-1) + numel(ptQuads_y(c,:));
    gd(counter_i:counter_f, c) = ptQuads_y(c, :);

end

% Embed the cut-out circle in the geometry
c = Nstrips+1;
gd(1, c) = 2;
gd(2, c) = numel(ptPoly_x);
gd(3:2+numel(ptPoly_x), c) = ptPoly_x(:);
gd(3+numel(ptPoly_y):end, c) = ptPoly_y(:);

%%
g = decsg(gd);
model = createpde;
geometryFromEdges(model,g);

%% Plot the geometry (for simple check)

figure(1);
pdegplot(model,'EdgeLabels','off')
axis equal


%% Mesh it
FEMesh = generateMesh(model,'Hmax', maxMeshSize, 'Hmin', ...
    minMeshSize, 'GeometricOrder', 'linear');
nodes2D = FEMesh.Nodes;

%% Identify which region is the main circular part
circlePartNo = 1;
minDist = outerRadius*10;
for c=1:Nstrips+1
    % Check if there is an element that is very close to zero
    cutOutElements = findElements(FEMesh,'region','Face', c);
    nodeNo = FEMesh.Elements(:, cutOutElements);
    xC = nodes2D(1, nodeNo);
    yC = nodes2D(2, nodeNo);
    dist = sqrt( xC.^2 + yC.^2 );
    if (minDist > min(dist(:)))
        circlePartNo = c;
        minDist = min(dist(:));
        % fprintf('Part number found: %d\n', c);
    end
end
cutOutElements = findElements(FEMesh,'region','Face', circlePartNo);


figure(2);
% F = pdegplot(model, 'FaceLabels','on');
pdeplot(model);
hold on
plot(nodes2D(1,:), nodes2D(2,:), 'ro');
cutOutElements = findElements(FEMesh,'region','Face',circlePartNo);
pdemesh(FEMesh.Nodes, FEMesh.Elements(:,cutOutElements),'EdgeColor','green');
hold off


return