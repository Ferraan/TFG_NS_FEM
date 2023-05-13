function Ne = StransfN (NeSCL, ndim)
%Maps from Scalar N to vector N
nnodeE = size (NeSCL, 2 ) ;
Ne = zeros ( ndim , ndim* nnodeE) ;
for inode = 1:nnodeE
    ini = ( inode-1)* ndim + 1; fin = inode * ndim ;
    Ne (:, ini : fin ) = NeSCL ( inode ) * eye ( ndim) ;
end
end