function Nst = AssemblyNGlobalP(Nelems,nelem,nnodeE,ndim,ngaus,CONNECT,nnode) ;
% This function returns the global "staked" N-matrix given the following
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
% J.A. Hernández, jhortega@cimne.upc.edu  (26 Oct. 2015, Barcelona, Spain)
if nargin == 0
    load('tmp1.mat')
end
m = nelem*ngaus ; % Number of rows 
n = nnode;%*ndim ;          % Number of columns
nzmaxLOC = m*nnodeE*ndim ;   % Maximum number of zeros (number of entries of Belems)
Nst = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Nst 
for inode = 1:nnodeE   % Loop over number of nodes at each element 
   setnodes = CONNECT(:,inode) ;  % Global numbering of the inode-th node of each element
   for idime = 1:1   % loop over spatial dimensions 
       % We employ here the Matlab's built-in "sparse" function --> S = sparse(i,j,s,m,n,nzmax). 
       %  This function uses
       % vectors i, j, and s to  generate an m-by-n sparse matrix
       % such that S(i(k),j(k)) = s(k), with space allocated
       % for nzmax nonzeros. 
       % --------------------------------------------------------
       % 1) Extracting the column of Nelems corresponding to the idime-th
       % DOF of the inode-th node of each element 
       DOFloc = (inode-1)*1+idime ;  % Corresponding indices
       s = Nelems(:,DOFloc) ;           % Corresponding column 
       % -----------------------------------------------------------------
       % 2) Now we determine the position within matrix Nst of each entry
       % of vector "s"  
       %  a) Rows
       i = [1:m]' ;  %    (all rows)
       %  b) Columns  
      DOFglo = (setnodes-1)+idime ;      %  Global DOFs associated to nodes   "setnodes" 
                                               %  and local DOF "idim". This vector has "nelem" entries,
                                               % while "i" has
                                               % "m=nstrain*nelem*ngaus ".
                                               % Therefore, we have to make
                                               % nstrain*ngaus tiling
                                               % copies of DOFglo (see next two lines)
       j = repmat(DOFglo', 1*ngaus,1) ; %     
       j = reshape(j,length(i),1) ;
        % Assembly process  (as sum of matrix with sparse representation)
       Nst = Nst + sparse(i,j,s,m,n,m) ; 
   end
end