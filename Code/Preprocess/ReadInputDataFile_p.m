function [COOR,CN,TypeElement,TypeElementB, ConductMglo,  rnod,dR,...  
    tracglo,CNb,fNOD,NameFileMesh]  = ReadInputDataFile_p(NameFileMeshDATA) 


% OUTPUTS 
% --------------
% 1. Finite element mesh 
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% TypeElementB: Type of boundary finite element (linear...)
% -----------
% 2. Material
% -----------
%  ConductMglo (ndim x ndim x nelem)  % Array of conductivity matrices
% -------------------------
% 3. Dirichlet (Essential) Boundary Condition s
% --------------------------------------------
%  rnod --> Set of nodes with prescribed temperatures 
%  dR   --> Vector of prescribed temperatures  (size(rnod) = size(dR))
% ---------------------------------------
% 4. Neumann (natural) Boundary Conditions
% -------------------------------------------
%  CNb: Connectivity matrix for the boundary elements of the Neumann
%  Boundary
%  qFLUXglo: Vector containing the prescribed flux at all nodes  (nnode x
% --------------------------------
% 5. Heat source
% ---------------
%  fNOD: Vector containing the nodal values of the heat source function (nnode x1 )
 disp('Reading input data...')
% Inputs example assigment 2 
% ----------------------------
%--------------------------------------------------------------------------
% 1. NAME OF THE MESH AND DATA FILES. COORDINATES, CONNECTIVITIES,  LISTS
% OF NODES FOR IMPOSING BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
% NameFileMeshDATA = 'Stokes30' ;   % For reading data from GID.

% -----------------------------------------------------------
% 2. Material data. 
% --------------------------------------------------------------
imat =1 ; % Index material
PROPMAT(imat).viscosity =  1  ; %  Kinematic viscosity "imat" (ISOTROPIC)
% -----------------------------------------------------------
% 3. Dirichlet boundary conditions (prescribed pressure)
% -----------------------------------------------------------

icond = 1; % Number of node condition 
DIRICHLET_Nodes(icond).NUMBER_NODE = 1 ;   % Number of node (for all nodes)
DIRICHLET_Nodes(icond).PRESCRIBED_P = 0 ;  % (constant along the line)


icond = 2; % Number of line condition 
DIRICHLET(icond).NUMBER_LINE = 1 ;   % Number of line
DIRICHLET(icond).PRESCRIBED_P = 0 ;  % (constant along the line)




% -------------------------------------------------
% 4. Neumann Boundary conditions (prescribed flux)
% ------------------------------------------------
icond= 1 ;
NEUMANN(icond).NUMBER_LINE = 1;  % Line 
NEUMANN(icond).PRESCRIBED_trac= 0 ;  % CONSTANT Prescribed traction


% -------------------------------------------
% 5. Momentum source (constant all over the body)
% --------------------------------------------
fSOURCE = 0 ; 


% END INPUTS 
% ----------------------------------------------------------------------------------------------







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB
NameFileMesh = [NameFileMeshDATA,'.msh']; % Name of the file containing the mesh information (Generated with GID)
% [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
%     ReadMeshFile(NameFileMesh)  ;
MESH=...
     ReadMeshFileStr_MULT(NameFileMesh,'READ_MATERIAL_COLUMN',1)  ;
 
 
 COOR = MESH.COOR; 
  CN = MESH.CN; 
  TypeElement = MESH.TypeElement;
  CONNECTb = MESH.CNb; 
  TypeElementB = MESH.TypeElementB; 
  MaterialType = MESH.MaterialType; 

nnode = size(COOR,1) ;% Number of nodes 

% 2. MATERIAL PROPERTIES: output ConductMglo  
%-----------------------
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements
ConductMglo = zeros(ndim,ndim,nelem) ; 
% Conductivity matrix (isotropic)
for imat = 1:length(PROPMAT)
    nu = PROPMAT(imat).viscosity ;
    ConductM = nu*eye(ndim) ; % eye = IDENTITY MATRIx
    ELEMS = find(MaterialType == imat) ;
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        ConductMglo(:,:,e) = ConductM ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which velocity is prescribed
%  Specify the number of line(s) in which the velocity is imposed
rnod = cell(length(DIRICHLET),1)  ; dR =  cell(length(DIRICHLET),1) ; 
for  icond = 1:length(DIRICHLET_Nodes)
    rnod{icond} =  ReadNodesMsh(NameFileMesh,DIRICHLET_Nodes(icond).NUMBER_NODE) ;
    dR{icond} = DIRICHLET_Nodes(icond).PRESCRIBED_P*ones(size( rnod{icond} )) ; 
end
for  icond = length(DIRICHLET_Nodes)+1:length(DIRICHLET)
    rnod{icond} =  ListOfNodesLINE(NameFileMesh,DIRICHLET(icond).NUMBER_LINE) ;
    dR{icond,1} = DIRICHLET(icond).PRESCRIBED_P*ones(size( rnod{icond} )) ; 
end
% Removed repeated condions 
% ---------------------------
rnod = cell2mat(rnod) ; 
dR = cell2mat(dR) ; 
[rnod, AAA] = unique(rnod) ;
dR = dR(AAA) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: qFLUXglo, CNb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CNb = cell(length(NEUMANN),1) ; % Boundary Connectivities 
tracglo =cell(length(NEUMANN),1) ;  % Prescribed, nodal values 

for icond = 1:length(NEUMANN)
    NODESb = ListOfNodesLINE(NameFileMesh,NEUMANN(icond).NUMBER_LINE) ;
    CNb{icond} = ElemBnd(CONNECTb,NODESb) ; 
    tracglo{icond} = NEUMANN(icond).PRESCRIBED_trac*ones(size(CNb{icond})) ;
end
 
CNb = cell2mat(CNb) ; 
tracglo = cell2mat(tracglo) ; 

%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Momentum source
%%
 
% Global vector of momentum sources (constant)
fNOD = fSOURCE*ones(nnode,1) ; 
