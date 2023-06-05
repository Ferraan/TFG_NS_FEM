function Ust = AssemblyUGlobal(u,nstrain,nelem,nnodeE,ndim,ngaus,CONNECT,nnode) ;
% This function returns the global "staked" U-matrix given the following
% inputs 
%  Nelems: N-matrices at all elements  (nelem*ngaus*nstrain x nnodeE*ndim)
% nstrain : Number of entries of the stress/strain vectors
% nelem: Number of elements
% nnodeE: Number of nodes per element 
% ndim: Number of spatial dimensions 
% ngaus: Number of gauss per element % 
% CONNECT: Array of connectivities 
% nnode: Number of nodes
% %

if nargin == 0
    load('tmp1.mat')
end
m = nstrain*nelem*ngaus ; % Number of rows 
n = nstrain*nelem*ngaus*2 ;          % Number of columns
nzmaxLOC = nelem*ngaus*4 ;   % Maximum number of zeros (number of entries of Belems)
Ust = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst 
for inode = 1:nnodeE   % Loop over number of nodes at each element 
   setnodes = CONNECT(:,inode) ;  % Global numbering of the inode-th node of each element
   for idime = 1:ndim   % loop over spatial dimensions 
       % % We employ here the Matlab's built-in "sparse" function --> S = sparse(i,j,s,m,n,nzmax). 
       % %  This function uses
       % % vectors i, j, and s to  generate an m-by-n sparse matrix
       % % such that S(i(k),j(k)) = s(k), with space allocated
       % % for nzmax nonzeros. 
       % % --------------------------------------------------------
       % % 1) Extracting the column of Belems corresponding to the idime-th
       % % DOF of the inode-th node of each element 
       % DOFloc = (inode-1)*ndim+idime ;  % Corresponding indices
       % s = Nelems(:,DOFloc) ;           % Corresponding column 
       % % -----------------------------------------------------------------
       % % 2) Now we determine the position within matrix Bst of each entry
       % % of vector "s"  
       % %  a) Rows
       % i = [1:m]' ;  %    (all rows)
       % %  b) Columns  
       DOFglo = (setnodes-1)*ndim+idime ;
       s=u(DOFglo); 
                                                %  Global DOFs associated to nodes   "setnodes" 
       %                                         %  and local DOF "idim". This vector has "nelem" entries,
       %                                         % while "i" has
       %                                         % "m=nstrain*nelem*ngaus ".
       %                                         % Therefore, we have to make
       %                                         % nstrain*ngaus tiling
       %                                         % copies of DOFglo (see next two lines)
       i=inode*nstrain-1:ngaus*nstrain:m;
       i2=i+1;
       j=(inode*nstrain*2-3+(-1+idime)*idime):nstrain*2*ngaus:n;
       j2=j+1;
       % j = repmat(DOFglo', nstrain*ngaus,1) ; %     
       % j = reshape(j,length(i),1) ;
       %  % Assembly process  (as sum of matrix with sparse representation)
       Ust = Ust + sparse(i,j,s,m,n,m) ;
       Ust = Ust + sparse(i2,j2,s,m,n,m);
   end
end
end