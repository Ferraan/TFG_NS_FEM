function Ahinv = MultiUpdateInverseHermitian(Binv,jrowMAT)
% Recursive application of UpdateInverseHermitian
%  J.A. Hdez, 24 March 2017
if nargin == 0
    m = 10; n =7; 
    jrowMAT =[8 9] ; 
    A = randn(m,n) ; a = randn(m,length(jrowMAT)) ;
    Bor = zeros(m,n+length(jrowMAT)) ; 
    indOR = 1:size(Bor,2) ; 
    indOR(jrowMAT) = [] ; 
    Bor(:,indOR) = A ; 
    Bor(:,jrowMAT) = a;     
    B = [Bor'*Bor] ;  
    Binv = inv(B);     
    AhinvREAL = inv(A'*A)      ;
        addpath('/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/FE_CODE/')

end
 
jrowMAT = sort(jrowMAT) ; 
BinvOLD = Binv ; 
for i = 1:length(jrowMAT)
    jrow = jrowMAT(i)-i+1 ; 
    Ahinv = UpdateInverseHermitian(BinvOLD,jrow) ; 
    BinvOLD = Ahinv ; 
end

