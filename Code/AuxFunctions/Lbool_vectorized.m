function Lbool = Lbool_vectorized (CN, nnode, ngaus , ndim)
nelem = size (CN, 1 ) ; nnodeE = size (CN, 2 ) ;
m = nelem* ndim* ngaus *nnodeE ; % Number of rows
n = nnode* ndim ; % Number of columns
nzmaxLOC = m*nnodeE* ndim ; % Maximum number of zeros
Lbool = sparse ( [ ] , [ ] , [ ] ,m, n , nzmaxLOC) ; % Allocating memory for Lbool
e = 1: nelem ;
s = ones ( size ( e ) ) ;
for igausLOC = 1: ngaus
    igausGLO = ( e-1)* ngaus + igausLOC ;
    for inodeLOC = 1: nnodeE
        inodeGLO = CN( e , inodeLOC ) ;
        irowGLO = ( igausGLO-1)* nnodeE + inodeLOC ;
        for idofLOC = 1: ndim
            idofGLO = ( irowGLO-1)* ndim + idofLOC ;
            for jdofLOC =1:ndim
                jdofGLO = ( inodeGLO-1)* ndim +jdofLOC ;
                if idofLOC == jdofLOC
                % Lbool_vector_gauss ( idofGLO , jdofGLO ) = 1;
                Lbool = Lbool + sparse ( idofGLO , jdofGLO , s ,m, n ,m) ;
                end
            end
        end
    end
end
end