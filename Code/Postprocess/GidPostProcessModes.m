function GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DOFl)

% Name of the mesh file 
NameFile_msh = [NameFileMesh(1:end),'_MODES.msh'] ; 
% Name of the results file 
NameFile_res= [NameFileMesh(1:end),'_MODES.res'] ; 

% Writing mesh file
MaterialType = ones(size(CN,1),1) ; 
GidMesh2DFE(NameFile_msh,COOR,CN,'MODES',MaterialType,TypeElement);
% Writing results file
%nMODES =size   ; 
MODESplot = zeros(size(COOR,1)*size(COOR,2),size(MODES,2)) ; 
MODESplot(DOFl,:) = MODES ; 
GidResults2DFE_modes(NameFile_res,COOR,CN,TypeElement,MODESplot,posgp);

cddd = cd ; 
NAMEFILEOPEN =  [cddd,'/',NameFile_res] ; 
disp('open GID FILE:')
disp(NAMEFILEOPEN)