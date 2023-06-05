function [C] = assemblyC(COOR_v,CN_v,u,TypeElement_v,Bst,OmegaGlo,debug,Nst)
%Assembly of matrix C
TypeIntegrand='K';
%disp('Begginning C assembly')

ndim = size(COOR_v,2);   % Spatial Dimension of the problem  (2 or 3)
nelem_v = size(CN_v,1);   % Number of elements 

nstrain=size(COOR_v,2); % Dimensions of the strain (velocity) != ndim always
ngaus=size(CN_v,2); %Number of Gauss points
nnode_v = size(COOR_v,1);  % Number of nodes
nnodeE_v = size(CN_v,2) ; %Number of nodes per element 



%%%%
time2=tic;
Ust=AssemblyUGlobal(u,nstrain,nelem_v,nnodeE_v,ndim,ngaus,CN_v,nnode_v);
%disp(['Time to assmble Ust fast: ',num2str(time1)]);
C=Nst'*Ust*OmegaGlo*Bst;
time2=toc(time2);
disp(['Time to assmble C fast: ',num2str(time2)]);
%%%%%%Naive method
if(debug)
    disp('Begginning C assembly slow method');
    [weig,~,shapef,dershapef] = ComputeElementShapeFun(TypeElement_v,nnodeE_v,TypeIntegrand);
    Cslow = sparse([],[],[],nnode_v*ndim,nnode_v*ndim,nnodeE_v*ndim*nelem_v) ;
    
    time1=tic;
    for e = 1:nelem_v
        %Nodes at element "e"
        CNloc_v = CN_v(e,:) ;          
        % Coordinates of the nodes of element "e"
        Xe_v = COOR_v(CNloc_v,:)' ;   
        
        % Computation of elemental C
        Ce = zeros(nnodeE_v*ndim,nnodeE_v*ndim);
        for  g = 1:length(weig)
            % Matrix of derivatives for Gauss point "g"
            BeXi_v = dershapef(:,:,g) ; 
            %Shapefunctions for pint "g" for scalars
            Nescalar_v=shapef(g,:); 
            %Velocity at node "e" and point "g"
            uel=u(2*CNloc_v(g)-1:2*CNloc_v(g),1);
            % Jacobian Matrix 
            Je_v = Xe_v*BeXi_v' ;
            
            % Jacobian 
            detJe = det(Je_v) ; 
            % Matrix of derivatives with respect to physical coordinates 
            BeTILDE_v = inv(Je_v)'*BeXi_v ;
            
            % Matrix of symmetric gradient 
            Be_v = QtransfB_special(BeTILDE_v,ndim);
            %Matrix Ne for vector functions
            Ne_v = StransfN(Nescalar_v,ndim);
            %Vector input
            u_vect=VtransfV(uel',ndim);
            % Element matrices
            Ce = Ce + weig(g)*detJe*(Ne_v'*u_vect*Be_v);
        end
  
        %Assembly of C
        for a = 1:nnodeE_v
           for b =1:nnodeE_v
               A = Nod2DOF(CN_v(e,a),ndim);
               B = Nod2DOF(CN_v(e,b),ndim);
               a1 = Nod2DOF(a,ndim);
               b1 = Nod2DOF(b,ndim);
               Cslow(A,B) = Cslow(A,B) + Ce(a1,b1);
           end
        end
        
        
        if (mod(e,1000) == 0)
           disp(strcat('e = ',int2str(e))) 
        end
    end
    time1=toc(time1);
    disp(['Time to assmble C: ',num2str(time1)]);
    if(~all(C-Cslow<1e-6,"all"))
                disp("Different")
    end
end
end
