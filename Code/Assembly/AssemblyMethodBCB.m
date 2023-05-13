function [K,Bst,OmegaGlo] = AssemblyMethodBCB(wSTs,Belem,nstrain,nelem,nnodeE,ndim,CN,nnode,wST,ViscMglo)
%%  Assembly  method in which the stiffness matrix is computed as: K = BstW^T*Cglo*BstW.
% ---------------------------------------------------------------------------
%dbstop('5')
    if nargin == 0
      load('tmp.mat')
    end
    
    ngaus = length(wSTs)/nelem ;
    disp('Assembly of Bst (stacked B-matrix)...')
    time2 = tic ; 
    Bst = AssemblyBGlobal(Belem,nstrain,nelem,nnodeE,ndim,ngaus,CN,nnode) ;
    % Bst is a ngaus*nelem*nstrain x ndofT matrix such that Est = Bst*d, where 
    % Est is a vector containing the strains at all Gauss points
    clear Belem
    time2 = toc(time2) ;
    disp(['Done (in =',num2str(time2),' s)']); 
    % Bst including weights
    % Diagonal matrix with weights
    %disp('Product BstW = W*Bst...')
    %mTOT = length(wST) ;
    %wDIAG = sparse(1:mTOT,1:mTOT,wST,mTOT,mTOT,mTOT) ;
    %BstW = wDIAG*Bst ;
    %disp('Done')
    % -----------------------
    % Matrix with all viscosity tensors (at all gauss points)
    time2 =tic ; 
    disp('Global viscosity matrix (including weights !!!)  ...')
    % Cglo =[C(e=1,g=1);C(e=1,g=2) ; ...
    Cglo_w = DefineElastMatGLO(ViscMglo,ngaus,wSTs)  ;
    % Cglo in block-diagonal format Cglo = diag(Cglo_1,Cglo_2 ...)
    Viscglo=ConvertBlockDiag(Cglo_w) ;
    OmegaGlo=ConvertBlockDiag(wST);
    time2 = toc(time2) ;
    disp(['Done (in =',num2str(time2),' s)']); 
    % % Stiffness matrix
    disp('Computing STIFFNESS MATRIX K = BstT.*Cglo_w*Bst  ...')
    %dbstop('34')
    time1 = tic ; 
    K = Viscglo*Bst ;
    K = (Bst)'*K;
    time1 = toc(time1) ; 
    disp(['Done (in =',num2str(time1),' s)']); 
end