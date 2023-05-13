function    detJ= determinantVECTORIZE(Jinp,ndim)
% Given a matrix Jinp = [Jinp_1; Jinp_22 ... ; Jinp_nelem]  (Jinp_i is ndim x ndim),
% determinantVECTORIZE returns a vector consisting of the determinants of
% each block matrix Jinp_i
% ndim = 2,3
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
if nargin == 0
    J1 = [3 5 6; 7 8 7 ; 0 0 3] ;  d1 =det(J1) ;
    J2 = 4*[3 5 6; 7 8 7 ; 0 0 3] ;  d2 =det(J2) ;
    J3 = 6*[3 5 6; 7 8 7 ; 0 0 3] ;  d3 =det(J3) ;
    Jinp = [J1; J2; J3];
    ndim = 3 ;
end
%

if   ndim == 1
    detJ = Jinp ;
else
    
    J = cell(ndim,ndim) ;
    for i=1:ndim
        iglo = i:ndim:size(Jinp,1) ;
        for j=1:ndim
            jglo = j ;
            J{i,j} = Jinp(iglo,jglo) ;
        end
    end
    if  ndim==2
        detJ = J{1,1}.*J{2,2}  -   J{1,2}.*J{2,1} ;
    elseif ndim == 3
        detJ = J{1,1}.*J{2,2}.*J{3,3} - J{1,1}.*J{2,3}.*J{3,2} - J{1,2}.*J{2,1}.*J{3,3} + J{1,2}.*J{2,3}.*J{3,1} ...
            + J{1,3}.*J{2,1}.*J{3,2} - J{1,3}.*J{2,2}.*J{3,1} ;
        
    end
    
end
