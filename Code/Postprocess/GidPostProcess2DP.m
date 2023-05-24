function GidPostProcess2DP(COOR,CN,TypeElement,p,NAME_INPUT_DATA,NameFileMesh) 
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
GidResults2DFEP(NameFile_res,COOR,CN,TypeElement,p,NAME_INPUT_DATA);

cddd = cd ; 
NAMEFILEOPEN =  [cddd,'/',NameFile_res] ; 
disp('open GID FILE FOR PRESSURE:')
disp(NAMEFILEOPEN)
end