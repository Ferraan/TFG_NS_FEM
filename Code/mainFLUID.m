%%% Finite Element Solver for the Navier-Stokes Flow  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Technical University of Catalonia
% Ferran de Miguel
% Conducted by Prof. Dr. Joaquin Hernández Ortega
% -------------------------------------------------------------------------


clear
close all

addpath("Preprocess/");
addpath("Assembly/");
addpath("Solver/")
addpath("Postprocess/");
addpath("AuxFunctions/");
%% INPUT  %% 
% Input data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAME_INPUT_DATA = 'DATA_ASSIGNMENT5' ;
%-------------------------------------------------


%% PREPROCESS 
% Inputs
% ---------
rho=1.225; %Air density [kg/m3]
mu=0.5; %Dynamic viscosty [kg/(m s)]
nu = mu/rho; %Kinematic Viscosity coefficient [m2/s] 
a0 = 0; %GLS stabilization parameter
Kovasznay=0; %Selects Boundary conditions for the kovasznay problem
NavierStokes=1; %Select 1 for NS, 0 for Stokes
res=1e-9; %Select residual for NS
maxiter=1e+3; %Maximum iterations for NS, if not converged exit with error
rel_factor=1; %Relaxation factor, <1 for under relaxation >1 for over relaxation
debug=0; %Debug parameter for enabling extra functions. 0 disabled, 1 enabled.

% ----------

% Definition of meshes 
% Name mesh
NameMeshP='LidDriven20'; %Cylinder10 LidDrivenBC5
NameMeshV=[NameMeshP '_v'];
% Velocity mesh
NameFileMeshDATA = ['./Meshes/' NameMeshV];
[COOR_v,CN_v,TypeElement_v,TypeElementB_v, ViscMglo,rnod_v,dR_v,...  
    ~,CNb_v,~,~] = ReadInputDataFile_v(NameFileMeshDATA,Kovasznay,nu); 
% Pressure mesh
NameFileMeshDATA = ['./Meshes/' NameMeshP];
[COOR_p,CN_p,TypeElement_p,TypeElementB_p, ~,rnod_p,dR_p,...  
    ~,CNb_p,~,~] = ReadInputDataFile_p(NameFileMeshDATA); 
% ----------

% Number of dimensions is obtained from the mesh


%% SOLVER
% Computation of K and G matrices
disp('Computing KGL')
[K,G,Bst,OmegaGlo] = ComputeKG(COOR_v,CN_v,COOR_p,CN_p,TypeElement_v, TypeElement_p,nu,debug,ViscMglo);

% Solution of the system of equations
if(NavierStokes)
    disp('Solving the system of equations for Navier Stokes')
    [u,v,p] = SolverNavierStokes(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,res,TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug);
else
    disp('Solving the system of equations for Stokes')
    [u,v,p] = SolverStokes(COOR_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G);
end

%% POST-PROCESS
disp('Starting the postprocessing')

%Postproc(COOR_v,COOR_p,u,v,p);
if(NavierStokes)
    NameMeshV=strcat(NameMeshV,'NavierStokesVisc_',num2str(nu));
    NameMeshP=strcat(NameMeshP,'NavierStokesVisc_',num2str(nu));
else
    NameMeshV=strcat(NameMeshV,'StokesVisc_',num2str(nu));
    NameMeshP=strcat(NameMeshP,'StokesVisc_',num2str(nu));
end
direc=['GIDPOST/',NameMeshP,'/'];
mkdir(direc);
GidPostProcess2DV(COOR_v,CN_v,TypeElement_v,u,v,NAME_INPUT_DATA,NameMeshV,direc); 
GidPostProcess2DP(COOR_p,CN_p,TypeElement_p,p,NAME_INPUT_DATA,NameMeshP,direc); 
