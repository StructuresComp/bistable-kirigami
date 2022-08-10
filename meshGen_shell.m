function meshGen_shell(nodesFile, bcFile,...
     stress, thickness_1, thickness_2,...
     FEMesh, cutOutElements, trilayer)

% This functions generates the 3D mesh of the structure from the 2D mesh input given

nodes2D = FEMesh.Nodes; % Nodes per sublayer
x2D = nodes2D(1, :);
y2D = nodes2D(2, :);
num_nodes2D = length(x2D);
el2D = FEMesh.Elements; % Elements per sublayer
[~,num_el2D] = size(el2D);

nLayers = 1; % Number of sub-layers per layer
zLayers = zeros(3*nLayers+1,1); % Z-coordinate of every sublayer
for i=1:3*nLayers+1
    if i==1
        zLayers(i)=0;
    elseif i<=nLayers+1
        zLayers(i)=(thickness_2/nLayers)*(i-1);
    elseif i<=2*nLayers+1
        zLayers(i) = thickness_2 + (thickness_1/nLayers)*(i-nLayers-1);
    elseif i<=3*nLayers+1
        zLayers(i) = thickness_2 + thickness_1 + (thickness_2/nLayers)*(i-2*nLayers-1);
    end
end

%%
fid = fopen(nodesFile, 'w'); %% write nodes and elements to file

%% Nodes
fprintf(fid, '*Node\n');

for cHeight=1:length(zLayers)
    for c2D = 1:num_nodes2D
        nodeNo = c2D + (cHeight-1)*num_nodes2D;
        z = zLayers(cHeight); % z = 0 for shell element

        fprintf(fid, '%d, %f, %f, %f\n', nodeNo, x2D(c2D), y2D(c2D), z);
    end
end
num_nodes = (3*nLayers+1)*num_nodes2D;
%% Elements
num_el_kir = numel(cutOutElements);

if trilayer == 1
    num_el = nLayers*(num_el2D + num_el_kir*2);
else
    num_el = nLayers*(num_el2D + num_el_kir);
end

el3D = zeros(6, num_el);

for cHeight = 1:nLayers    % Bottom Kirigami Layer
    for c2D = 1:numel(cutOutElements)

        cEl = cutOutElements(c2D);% element number

        elNo = c2D+num_el_kir*(cHeight-1);

        offset = num_nodes2D*(cHeight-1);

        el3D(1, elNo) = el2D(1, cEl) + offset;
        el3D(2, elNo) = el2D(2, cEl) + offset;
        el3D(3, elNo) = el2D(3, cEl) + offset;
        el3D(4, elNo) = el2D(1, cEl) + offset + num_nodes2D;
        el3D(5, elNo) = el2D(2, cEl) + offset + num_nodes2D;
        el3D(6, elNo) = el2D(3, cEl) + offset + num_nodes2D;
    end
end

for cHeight = 1:nLayers   %Substrate Layer
    for c2D = 1:num_el2D
        elNo = c2D +num_el2D*(cHeight-1)+nLayers*num_el_kir;

        offset = num_nodes2D*nLayers+num_nodes2D*(cHeight-1);        

        el3D(1, elNo) = el2D(1, c2D) + offset;
        el3D(2, elNo) = el2D(2, c2D) + offset;
        el3D(3, elNo) = el2D(3, c2D) + offset;
        el3D(4, elNo) = el2D(1, c2D) + offset + num_nodes2D;
        el3D(5, elNo) = el2D(2, c2D) + offset + num_nodes2D;
        el3D(6, elNo) = el2D(3, c2D) + offset + num_nodes2D;        
    end
end

if trilayer == 1 % Top Kirigami Layer (only in case of trilayer)
    for cHeight = 1:nLayers    
        for c2D = 1:numel(cutOutElements)

            cEl = cutOutElements(c2D); % element number

            elNo = c2D + num_el_kir*(cHeight-1)+nLayers*num_el_kir+nLayers*num_el2D;

            offset = num_nodes2D*2*nLayers+num_nodes2D*(cHeight-1);

            el3D(1, elNo) = el2D(1, cEl) + offset;
            el3D(2, elNo) = el2D(2, cEl) + offset;
            el3D(3, elNo) = el2D(3, cEl) + offset;
            el3D(4, elNo) = el2D(1, cEl) + offset + num_nodes2D;
            el3D(5, elNo) = el2D(2, cEl) + offset + num_nodes2D;
            el3D(6, elNo) = el2D(3, cEl) + offset + num_nodes2D;
        end
    end
end

fprintf(fid, '*Element, type=C3D6\n'); %% 3-node triangular shell element
for c = 1:num_el
        fprintf(fid, '%d, %d, %d, %d, %d, %d, %d \n', c, el3D(1, c),el3D(2, c), ...
            el3D(3, c),el3D(4, c),el3D(5, c), ...
            el3D(6, c));
end

%% Element and Node sets

fprintf(fid,'*Nset, nset = NodeBot, generate\n'); % all nodes in the bottom layer
fprintf(fid,'%d, %d, 1\n', 1, (nLayers+1)*num_nodes2D);
fprintf(fid, '*Elset, elset=elSubstrate, generate\n'); % all elements in the substrate layer
fprintf(fid, '%d, %d, 1\n', num_el_kir*nLayers+1, num_el_kir*nLayers+num_el2D*nLayers);

fprintf(fid,'*Nset, nset = NodeSubstrate, generate\n'); % all nodes in the substrate layer
fprintf(fid,'%d, %d, 1\n', nLayers*num_nodes2D+1, (2*nLayers+1)*num_nodes2D);
fprintf(fid, '*Elset, elset = elBot, generate\n'); % all elements in the substrate layer

fprintf(fid,'%d, %d, 1\n',1, num_el_kir*nLayers);
if trilayer==1
    fprintf(fid,'*Nset, nset = NodeTop, generate\n'); % all nodes in the top layer
    fprintf(fid,'%d, %d, 1\n', 2*nLayers*num_nodes2D+1, num_nodes);
    fprintf(fid, '*Elset, elset = elTop, generate\n'); % all elements in the top layer
    fprintf(fid,'%d, %d, 1\n',num_el_kir*nLayers+num_el2D*nLayers+1, num_el);
end

fprintf(fid,'*Nset, nset = NodeAll, generate\n'); % all nodes
fprintf(fid,'%d, %d, 1\n', 1, num_nodes);

%  Center Node
% Find node located at (0,0) at the center of the composite
xDesired = 0.0;
yDesired = 0.0;
dist = sqrt( (x2D - xDesired).^2 + (y2D - yDesired).^2 );
[~, minInd] = min(dist);
fprintf(fid, '*Nset, nset=centerNode\n%d\n', minInd+num_nodes2D*(floor((nLayers*3+1)/2-1)));

fprintf(fid, '*Solid Section, elset=elSubstrate, material=Material-Substrate\n');
fprintf(fid, '*Solid Section, elset=elBot, material=Material-Kirigami\n');
fprintf(fid, '*Solid Section, elset=elTop, material=Material-Kirigami\n');     

fclose(fid);


%% Boundary Conditions and Stresses

fid = fopen(bcFile,'w');
% Fix the center node of the composite
fprintf(fid,'**Name: BC1 Type: Displacement/Rotation\n');
fprintf(fid,'*Boundary\n');
fprintf(fid, 'centerNode, 1, 6,%f\n', 0.0);

% Initiate Pre-stretch in terms of initial stress conditions
fprintf(fid,'**\n** PREDIFINED FIELDS\n**\n');
fprintf(fid,'*Initial Conditions, type = STRESS\n');
fprintf(fid,'elSubstrate, %f, %f, %f\n', stress, stress, 0.0); 

fclose(fid);
