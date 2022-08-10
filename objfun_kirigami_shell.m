function objfun_kirigami_shell(main_dimension, lambda, stress, ...
     thickness_1, thickness_2,mat_param,...
     FEMesh, cutOutElements,...
     abaqus_addr, diagnostic, trilayer)

%%
critTime = Inf; % If the Abaqus simulation does not finish within
% "critTime" (seconds), it is forcefully closed

checkInt = 30; % Check status of Abaqus simulation every few seconds 
% specified by this variable.


%% STEP 1: Generate mesh

if trilayer==1
    inputNameString = 'r-tri-';
else
    inputNameString='r-bi-';
end

% Name of the job
jobName = sprintf( [inputNameString, '%d'], main_dimension*1000);

% Name of the input file. It uses template.inp (or some other "base"
% input file as the foundation and then copies+edits it. 
% See write_input_file() function for details.
inputName = sprintf( [inputNameString, '%d-%d.inp'], main_dimension*1000, lambda*100);

% The input file requies TWO more files that will need to be written by
% MATLAB:
% (1) nodesFile (nodes and elements of substrate layer)
% (2) bcfile (contains boundary conditions)

% Name of the file containing the nodes and elements of the substrate layer
nodesFile = sprintf( 'r-nodes-%d-%d.inp', main_dimension*1000, lambda*100);

% Boundary condition file
bcFile = sprintf( 'r-boundaryCondition-%d-%d.inp', main_dimension*1000, lambda*100);
%% Create an input (.inp) file
write_input_file(inputName, nodesFile, bcFile, mat_param); 

%% Populate the files required
meshGen_shell(nodesFile, bcFile,...
    stress,thickness_1, thickness_2,...
     FEMesh, cutOutElements, trilayer);
%%
Cmd1 = [abaqus_addr,' j=', jobName, ' input=', inputName, ' &'];
Cmd2 = [abaqus_addr,' terminate job=', jobName];

%% STEP 2: Run simulation
system(Cmd1);
pause(5);
%%
% Something to keep in mind:
% For some reason, the abaqus software on the desktop does not close by
% itself. The remaining of this code section deals with this problem. It
% checks the status using .sta file every few seconds (variable checkInt).
% If this was not problem, we would simply need the following two commands
% instead of this entire section.
% Cmd = [abaqus_addr,' j=', jobName, ' input=', inputName];
% system(Cmd);


% Check for *sta file
staFile = [jobName, '.sta'];

jobComplete = false;
odbname = [jobName, '.odb'];

ind = [];

runTime = 0;
while (jobComplete == false)
    
    % If file doesn't exist, let us continue
    if isfile(staFile) == false
        continue
    end
    
    fid = fopen(staFile, 'r');
    while ~feof(fid)
        tline = fgetl(fid);
    end
    fclose(fid);
    
    ind = strfind(tline, 'HAS COMPLETED SUCCESSFULLY');
    ind2 = strfind(tline, 'HAS NOT BEEN COMPLETED');

    if numel(ind)>0
       jobComplete = true;
       fprintf('THE ANALYSIS HAS COMPLETED SUCCESSFULLY\n');
        system(Cmd2);
    %   system(['fuser -k ', jobName, '*']); %If running multiple jobs and want to close all processes to save space

    elseif numel(ind2)>0
        jobComplete = true;
        fprintf('THE ANALYSIS HAS NOT BEEN COMPLETED\n');
        system(Cmd2);
    %   system(['fuser -k ', jobName, '*']); %If running multiple jobs and want to close all processes to save space

    elseif runTime > critTime
        jobComplete = true;
        fprintf('THE ANALYSIS HAS NOT BEEN COMPLETED AS THE SIMULATION TOOK TOO LONG\n');
        system(Cmd2);
    %   system(['fuser -k ', jobName, '*']); %If running multiple jobs and want to close all processes to save space

    else
        jobComplete = false;
    end

        pause(checkInt);
        runTime = runTime + checkInt;
end
%
% At this point, if "ind" is empty, the simulation did not finish
% successfully. Otherwise, if numel(ind)>0, the simulation finished
% successfully.
%

% Command to read the odb to get required outputs

Cmd3 = [abaqus_addr,' cae noGUI=readODB -- ', odbname];

if numel(ind)>0 % successful completion
    %% STEP 3: Read ODB file
    if isfile(odbname) == true
        system(Cmd3);
    else
        fprintf('%s file does not exist\n', odbname);
    end
    
end

%% Delete files
switch diagnostic
    case 1              
        system('rm abaqus*');
        system(['rm ', jobName, '*']);
        system(['rm ', inputName, ' ', nodesFile, ' ', bcFile]);
    case 2
        system(['mv ', odbname, ' final-', odbname])
        system('rm abaqus*');
        system(['rm ', jobName, '*']);
        system(['rm ', inputName, ' ', nodesFile,' ', bcFile]);
    otherwise
       
end

end
