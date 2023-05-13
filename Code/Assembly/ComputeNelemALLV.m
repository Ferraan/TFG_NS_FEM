function[ Nelem] = ComputeNelemALLV(COOR,CN,TypeElement,nstrain) 

    nnode = size(COOR,1); ndim =2; nelem = size(CN,1); nnodeE = size(CN,2) ;
    % nstrain = size(celasglo,1) ;
    % Shape function routines (for calculating shape functions and derivatives)
    TypeIntegrand = 'K'; %Although its C, the matrices are the same
    %Calculate element weights and shape functions
    [weig,~,shapef,~] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
    
    
    ngaus = length(weig) ;
    % Initialization
    
    Nelem = zeros(nstrain*nelem*ngaus,nnodeE*2) ;
    % COORDINATE MATRIX arranged in a nelem*ndim x nnodeE matrix
    % Let us define a matrix ROWSgauss such that   Nelem(ROWSgauss(g,:),:) returns 
    % the N-matrices of the g-th points of all elements
    indREF = 1:nstrain*ngaus*nelem ;
    ROWSgauss = reshape(indREF,nstrain,nelem*ngaus) ; 
    for  g = 1:ngaus
        
        
        ROWSglo =  ROWSgauss(:,g:ngaus:ngaus*nelem);   
        ROWSglo = ROWSglo(:) ;
       
        NeALL=repmat(shapef(g,:),nelem,1);
        NeALL=QtransfBvect(NeALL,2);
        Nelem(ROWSglo,:)=NeALL;
         
        
    end

end