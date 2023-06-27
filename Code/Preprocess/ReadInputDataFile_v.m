function [COOR,CN,TypeElement,TypeElementB, ViscMglo,  rnod,dR,...  
    qFLUXglo,CNb,fNOD,NameFileMesh]  = ReadInputDataFile_v(NameFileMeshDATA,Ux,nu) 

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
%  rnod --> Set of nodes with prescribed velocities
%  dR   --> Vector of prescribed velocities  (size(rnod) = size(dR))
% ---------------------------------------
% 4. Neumann (natural) Boundary Conditions
% -------------------------------------------
%  CNb: Connectivity matrix for the boundary elements of the Neumann
%  Boundary
% --------------------------------
% 5. Heat source
% ---------------
%  fNOD: Vector containing the nodal values of the traction function (nnode x1 )
disp('Reading input data for velocity...')

% Inputs example assigment 2 
% ----------------------------
%--------------------------------------------------------------------------
% 1. NAME OF THE MESH AND DATA FILES. COORDINATES, CONNECTIVITIES,  LISTS
% OF NODES FOR IMPOSING BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
% -----------------------------------------------------------
% 2. Material data. 
% --------------------------------------------------------------
imat =1 ; % Index material
PROPMAT(imat).viscosity =  nu  ; %  Kinematic viscosity "imat" (ISOTROPIC)
% -----------------------------------------------------------
% 3. Dirichlet boundary conditions (prescribed velocity)
% -----------------------------------------------------------
icond = 1; % Number of condition 
DIRICHLET(icond).NUMBER_LINE = 1 ;   % Number of line
DIRICHLET(icond).PRESCRIBED_Ux = 0 ;  % (constant along the line)
DIRICHLET(icond).PRESCRIBED_Uy = 0 ;
icond = 2; % Number of condition 
DIRICHLET(icond).NUMBER_LINE = 2 ;   
DIRICHLET(icond).PRESCRIBED_Ux = Ux ;
DIRICHLET(icond).PRESCRIBED_Uy = 0 ;



% -------------------------------------------------
% 4. Neumann Boundary conditions (prescribed flux)
% ------------------------------------------------
icond= 1 ;
NEUMANN(icond).NUMBER_LINE = 1;  % Line 
NEUMANN(icond).PRESCRIBED_qBAR= 0 ;  % CONSTANT Prescribed heat flux vector x   normal unit vector to the line


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
ViscMglo = zeros(ndim*2,ndim*2,nelem) ; 
% Conductivity matrix (isotropic)
for imat = 1:length(PROPMAT)
    nu = PROPMAT(imat).viscosity ;
    ConductM = nu*eye(ndim*2) ; % eye = IDENTITY MATRIx
    ELEMS = find(MaterialType == imat) ;
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        ViscMglo(:,:,e) = ConductM ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which velocity is prescribed
%  Specify the number of line(s) in which the velocity is imposed
rnod = cell(length(DIRICHLET),1)  ; dR =  cell(length(DIRICHLET),1) ; 
for  icond = 1:length(DIRICHLET)
    rnod{icond} =  ListOfNodesLINE(NameFileMesh,DIRICHLET(icond).NUMBER_LINE) ;
    dR{icond,1} = DIRICHLET(icond).PRESCRIBED_Ux*ones(size( rnod{icond} )) ; 
    dR{icond,2} = DIRICHLET(icond).PRESCRIBED_Uy*ones(size( rnod{icond} )) ; 
end
% Removed repeated condions 
% ---------------------------
rnod = cell2mat(rnod) ; 
rnod=[2*rnod-1; 2*rnod];
dR = cell2mat(dR) ; 
[rnod, AAA] = unique(rnod) ;
dR = dR(AAA) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: qFLUXglo, CNb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CNb = cell(length(NEUMANN),1) ; % Boundary Connectivities 
qFLUXglo =cell(length(NEUMANN),1) ;  % Prescribed, nodal values 

for icond = 1:length(NEUMANN)
    NODESb = ListOfNodesLINE(NameFileMesh,NEUMANN(icond).NUMBER_LINE) ;
    CNb{icond} = ElemBnd(CONNECTb,NODESb) ; 
    qFLUXglo{icond} = NEUMANN(icond).PRESCRIBED_qBAR*ones(size(CNb{icond})) ;
end
 
CNb = cell2mat(CNb) ; 
qFLUXglo = cell2mat(qFLUXglo) ; 

%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Momentum source
%%
 
% Global vector of momentum sources (constant)
fNOD = fSOURCE*ones(nnode,1) ; 
