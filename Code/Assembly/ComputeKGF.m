function [K,G,F,Bst,OmegaGlo] = ComputeKGF(COOR_v,CN_v,COOR_p,CN_p,CNb,tracglo,TypeElement_v, TypeElement_p,TypeElementB,nu,ViscMglo,debug) 
%%%%
% This subroutine   returns the global conductance matrix K,
% the global gradient operator G
% B stacked matrix Bst for velocity nodes 
% OmegaGlo matrix containing the global weights at each Gauss point
% Inputs
% --------------
% 1. Finite element mesh
% -------------------
% COOR_v: Coordinate matrix of velocities (nnode x ndim) 
% CN_v:   Connectivity matrix of veloicities (nelem x nnodeE)
% TypeElement_v: Type of finite element for velocities(quadrilateral,...)
% COOR_p: Coordinate matrix of pressure (nnode x ndim) 
% CN_p:  Connectivity matrix of pressure (nelem x nnodeE)
% TypeElement_p: Type of finite element for pressure (quadrilateral,...)
% nu: fluid kinematic viscosty 
% debug: parameter for controlling the assembly process
% ViscMglo: global viscosity matrix containing the viscosity at each
% element 
%%%%


 
% Dimensions of the problem
ndim = size(COOR_v,2);   % Spatial Dimension of the problem  (2 or 3)
nelem = size(CN_v,1);   % Number of elements 

nstrain=size(COOR_v,2); % Dimensions of the strain (velocity) != ndim always

nnode_v = size(COOR_v,1);  % Number of nodes
nnodeE_v = size(CN_v,2) ; %Number of nodes per element 

nnode_p = size(COOR_p,1);  % Number of nodes
nnodeE_p = size(CN_p,2) ; %Number of nodes per element 



%%%Assembly of K
[ Belem, wSTs, wST, ~] = ComputeBelemALL(COOR_v,CN_v,TypeElement_v,nstrain) ;
[K,Bst,OmegaGlo] = AssemblyMethodBCB(wSTs,Belem,nstrain,nelem,nnodeE_v,ndim,CN_v,nnode_v,wST,ViscMglo) ;
%%%Assembly of G

[Nelem]=ComputeNelemALL(COOR_p,CN_p,TypeElement_p,1);
[G]=AssemblyNIB(Nelem,nelem,4,ndim,9,CN_p,nnode_p,OmegaGlo,Bst);

% Assembly of global traction vector (boundary contribution)
disp('Computing  global tractions vector (boundary contribution) Fbnd...')
Fbnd = FdisCOMP(COOR_v,CNb,TypeElementB, tracglo,2) ; 
%%%%Body forces not coded.
F=Fbnd;

%%%%Naive method
    if(debug)
        disp('Begginning K and G naive method assembly');
        % Determine Gauss weights, shape functions and derivatives  
        TypeIntegrand = 'K'; 
        time1=tic;
        Kslow = sparse([],[],[],nnode_v*ndim,nnode_v*ndim,nnodeE_v*ndim*nelem) ;
        
        Gslow = sparse([],[],[],nnode_p,nnode_v*ndim,nnodeE_v*ndim*nelem) ;
        
        [weig_v,~,shapef_v,dershapef_v] = ComputeElementShapeFun(TypeElement_v,nnodeE_v,TypeIntegrand) ; 
        [~,~,shapef_p,dershapef_p] = ComputeElementShapeFun(TypeElement_p,nnodeE_p,TypeIntegrand) ; 
        for e = 1:nelem
            
            CNloc_v = CN_v(e,:) ;          % Coordinates of the nodes of element "e"
            Xe_v = COOR_v(CNloc_v,:)' ;      % Computation of elemental stiffness matrix
            
            CNloc_p = CN_p(e,:) ;          % Coordinates of the nodes of element "e"
            Xe_p = COOR_p(CNloc_p,:)' ;      % Computation of elemental stiffness matrix
            
            
            % Compute element matrix
            % ----------
            [Ke,Ge,~] = ComputeKeGeLeMatrix(weig_v,shapef_v,dershapef_v,shapef_p,dershapef_p,Xe_v,Xe_p,nu) ;
            % ----------
            % Assembly K
            for a = 1:nnodeE_v
               for b =1:nnodeE_v
                   A = Nod2DOF(CN_v(e,a),ndim);
                   B = Nod2DOF(CN_v(e,b),ndim);
                   a1 = Nod2DOF(a,ndim);
                   b1 = Nod2DOF(b,ndim);
                   Kslow(A,B) = Kslow(A,B) + Ke(a1,b1);
               end
            end
            % Assembly G
            for a = 1:nnodeE_v
               for b =1:nnodeE_p
                   A = Nod2DOF(CN_v(e,a),ndim);
                   B = CN_p(e,b);
                   a1 = Nod2DOF(a,ndim);
                   Gslow(B,A) = Gslow(B,A) + Ge(b,a1);
               end
            end    
            
            
            
            if (mod(e,1000) == 0)
               disp(strcat('e = ',int2str(e))) 
            end
            end
            %By definition
            Gslow = -Gslow;
            %%validation
            time1=toc(time1);
            disp(['Time to assmble K and G slow: ',num2str(time1)]);
            if(~all(Kslow-K<1e-6,"all") || ~all(Gslow-G<1e-6,"all"))
              warning("Different")
            end
        
    end
end

