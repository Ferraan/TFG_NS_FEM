function[ Nelem] = ComputeNelemALL(COOR,CN,TypeElement,nstrain) 
%%%%
% This subroutine   returns
%  1) the matrix Nelem (nstrain*nelem*ngaus x ndime*nnodeE)  consisting of
%  all element N-matrices for pressure nodes
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE),
% TypeElement: Type of finite element (quadrilateral,...), % nstrain dimensions of the strain 
%% Vectorized version
nnode = size(COOR,1); ndim = 1; nelem = size(CN,1); nnodeE = size(CN,2) ;
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'G';
%dbstop('20')
[weig,~,shapef,~] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
ngaus = length(weig) ;
% Initialization
Nelem = zeros(nstrain*nelem*ngaus,nnodeE) ;
% COORDINATE MATRIX arranged in a nelem*ndim x nnodeE matrix
% Let us define a matrix ROWSgauss such that   Belem(ROWSgauss(g,:),:) returns 
% the N-matrices of the g-th points of all elements
indREF = 1:nstrain*ngaus*nelem ;
ROWSgauss = reshape(indREF,nstrain,nelem*ngaus) ; 
for  g = 1:ngaus
    
    
    ROWSglo =  ROWSgauss(:,g:ngaus:ngaus*nelem);   
    ROWSglo = ROWSglo(:) ;
   
    NeALL=repmat(shapef(g,:),nelem,1);
    Nelem(ROWSglo,:)=NeALL;
     
    
end

end