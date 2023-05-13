function [u_vect] = VtransfV(u,ndim)
%Transfer function from 1x2 to 2x4
nnodeE = size (u, 2 ) ;
u_vect = zeros ( ndim , ndim* nnodeE) ;
for inode = 1:nnodeE
    ini = ( inode-1)* ndim + 1; fin = inode * ndim ;
    u_vect (:, ini : fin ) = u ( inode ) * eye ( ndim) ;
end
end