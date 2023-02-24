%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    This code is for soft kirigami simulation  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

trilayer = 1; % 1 for trilayer, else for bilayer
abaqus_addr = '/home/sci02/abaqus/Commands/abaqus'; % Location to abaqus on the computer

%% Kirigami Shape
shape_id = 1; % 1 for Cross Shape, 2 for Lotus, you can add your shape and mesh file here.

switch shape_id
    case 1
        % Inputs parameters for Cross Shape
        L = 40/1000;  % Length of the Cross Shape (m)
        H=0.005*L/0.06; % Height of the Cross Shape (m)
        Xk=H/2; % X coordinate required for mesh definition
        Yk=sqrt((L/2)^2-Xk^2); % Y coordinate
        main_dimension = L/2;
        maxMeshSize = main_dimension/20; % max size of mesh
        minMeshSize = maxMeshSize/2; % min size of mesh
        
        [FEMesh, cutOutElements] = cutoutMesh_cross(L, Xk, Yk,maxMeshSize,minMeshSize);
    case 2
        % Inputs parameters for Lotus Shape
        R = 20/1000;  % OuterRadius (m)
        radiusRatio = 0.550; % Ratio of the inner Radius to the Outer Radius
        Nstrips = 10; % Number of lotus petal-lls
        innerRadius = radiusRatio*R;        
        alpha = 0.4; % Percentage Material removed
        dTheta = 2*pi*alpha/Nstrips; % Angle required for mesh definition - central angle for cutout
        
        main_dimension = R;
        maxMeshSize = main_dimension/20; % max size of mesh
        minMeshSize = maxMeshSize/2; % min size of mesh
        
        [FEMesh, cutOutElements] = cutoutMesh_lotus(innerRadius,R,...
            Nstrips, dTheta, maxMeshSize, minMeshSize);

    otherwise 
        disp("Please define shape")
        exit()        
end

diagnostic = 0; % 1 to delete all abaqus files post simulation, 2 to delete all except odb, otherwise not delete any files
%%
lambda = 1.2; % Substrate prestretch
thickness_1 = (1.1/lambda^2) / 1000; % Thickness of the substrate | reduce thickness b/c of post-stretching, divided by 2-2*lambda
thickness_2 = 1.6 / 1000; % Thickness of the kirigami

%% material properties
% Mooney-Rivlin Parameters
% Substrate
C1s = 2.206e4; 
C2s = 1656;
D1s = 0.00001;

%Kirigami
C1k = 1.788e4;
C2k = 8.445e4;
D1k = 0.00001;

mat_param = [C1s, C2s, D1s, C1k, C2k, D1k];

%% Stress Calculation
I1 = 2*lambda^2 + 1/lambda^4; % stress invariant
dI1 = 4*lambda - 4/lambda^5; % derivative of stress invariant
stress = 2*C1s*(lambda^2-1/lambda^4)-2*C2s*(1/lambda^2-lambda^4); % stress calculation from pre-stretch

%% Build and start the simulation

objfun_kirigami_shell(main_dimension, lambda, stress, ...
    thickness_1, thickness_2, mat_param,...
    FEMesh, cutOutElements,...
    abaqus_addr, diagnostic, trilayer);



