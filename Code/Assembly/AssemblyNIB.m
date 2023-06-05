function [G]=AssemblyNIB(Nelems,nelem,nnodeE,ndim,ngaus,CONNECT,nnode_p,OmegaGlo,Bst) 
%Computes G as the product of G = NstW^T*F*OmegaGlo*BstW.
%Inputs 
   % Nelems=matrix of pressure shapefunctions of dims ngaus*nelems x 4
   %    nstrain= global strain (2) for vx and vy
   %    nelem= total number of elements
   %    nnodeE= number of pressure nodes per element
   %    ndim= number of dimensions
   %    ngaus= number of gauss points
   %    CONNECT= connectivity matrix for pressure elements
   %    nnode_p= number of total pressure nodes
   %    OmegaGlo= matrix containing the product of all gauss weights and
   %    jacobians
   %    Bst= B stacked matrix for velocity
%Compute the N stacked matrix
 time2 =tic ; 
NstP=AssemblyNGlobalP(Nelems,nelem,nnodeE,ndim,ngaus,CONNECT,nnode_p);
%Compute matrix F as the Kronecker tensor product
f=speye(nelem*ngaus);
I=kron(f,[1 0 0 1]);
%Assembly of G
G=NstP'*I*OmegaGlo*Bst;
G=-G;
time2 = toc(time2) ;
disp(['Time to assemble G fast ',num2str(time2),' s']); 
end