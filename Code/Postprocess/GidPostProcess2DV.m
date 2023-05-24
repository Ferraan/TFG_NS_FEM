function GidPostProcess2DV(COOR,CN,TypeElement,u,v,NAME_INPUT_DATA,NameFileMesh) 
% Post-processing of results using GID (2D)
if nargin==0
    load('tmp.mat')
end
% Name of the mesh file 
NameFile_msh = [NameFileMesh,'.msh'] ; 
% Name of the results file 
NameFile_res= [NameFileMesh,'.res'] ; 

% Writing mesh file
MaterialType = ones(size(CN,1),1) ; 
GidMesh2DFE(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement);
% Writing results file
GidResults2DFEV(NameFile_res,COOR,CN,TypeElement,u,v,NAME_INPUT_DATA);

cddd = cd ; 
NAMEFILEOPEN =  [cddd,'/',NameFile_res] ; 
disp('open GID FILE FOR VELOCITY:')
disp(NAMEFILEOPEN)
end