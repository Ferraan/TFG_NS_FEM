function [u,v,p] = SolverNavierStokes(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,res,TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug)
%Inputs:
%   - COOR_v: Matrix with the coordiantes of velocity nodes
%   - COOR_p: Matrix with the coordiantes of pressure nodes
%   - rnod_v: Vector containing the restricted nodes of the velocity
%   - rnod_p: Vector containing the restricted nodes of the pressure
%   - dR_v: Vector containing the prescribed velocities
%   - dR_p: Vector containing the prescribed pressures
%   - K: Global viscosity matrix 
%   - G: Global gradient operator
%   - res: tolerance of the residual for the nonlinear system of equations
%   - TypeElement_v: Velocity type of element for computing C
%   - maxiter: maximum iterations before the nonlinear solver returns an
%   error
%   - rel_factor: relaxation factor for the nonlinear solver
%   - Bst: B-stacked matrix for velocity
%   - OmegaGlo: Matrix containing the products of all Jacobians and weights
%   at all Gauss points
%   - debug: Debug parameter to control the assembly of C
%Outputs:
%   - u: vector with the horizontal component of the velocity
%   - v: vector with the vertical component of the velocity
%   - p: vector with the pressure
    ndim = size(COOR_v,2); 
  
    nDOF_v = ndim * size(COOR_v,1);
    DOFl_v = (1:nDOF_v)';
    DOFr_v = rnod_v;
    DOFl_v(DOFr_v) = [];
    
    % Pressure on node 1 set to zero
    nDOF_p = size(COOR_p,1);
    DOFl_p = (1:nDOF_p)';
    DOFr_p = rnod_p;
    DOFl_p(DOFr_p) = [];
    
    
    % Decomposition of matrices
    Kll = K(DOFl_v,DOFl_v);
    Gl = G(DOFl_p,DOFl_v);
    GlT = Gl';
    Gr = G(DOFl_p,DOFr_v);
    Klr = K(DOFl_v,DOFr_v);
    %Auxiliary matrix
    L = zeros(length(DOFl_p));

   
    %Precalculation of matrix Nst for velocity
    nnode_v = size(COOR_v,1); ndim = size(COOR_v,2); nelem_v = size(CN_v,1); nnodeE_v = size(CN_v,2) ; ngaus=size(CN_v,2); nstrain=size(COOR_v,2);
    Nelem=ComputeNelemALLV(COOR_v,CN_v,TypeElement_v,2);
    Nst=AssemblyNGlobalV(Nelem,nstrain,nelem_v,nnodeE_v,ndim,ngaus,CN_v,nnode_v);
    
    %Solve of the nonlinear system
    disp('Solving Navier-Stokes...')
    % Allocate ones for first iteration
    u=ones(nDOF_v,1);
    u(DOFr_v)=dR_v;
    u_new=zeros(nDOF_v,1);
    normResidual=1; % Inicialization of the residual
    iter=0;
    time1=tic;
    while (res<normResidual && iter<maxiter)
        %Assembly of matrix C, recall that C(u)
        C=assemblyC(COOR_v,CN_v,u,TypeElement_v,Bst,OmegaGlo,debug,Nst);
        Cll=C(DOFl_v,DOFl_v);
        Clr=C(DOFl_v,DOFr_v);
        %Solve the linear system
        d = [Kll+Cll GlT; Gl L]\[-Klr*dR_v+Clr*dR_v; -Gr*dR_v];
        %Assign the corresponding DOF to the velocity
        u_new(DOFl_v)=d(1:length(DOFl_v));
        u_new(DOFr_v)=dR_v;
        %Calculate the residual
        normResidual=max(norm(u_new-u));
        %Update the velocity with a relaxation factor
        u=(1-rel_factor)*u+rel_factor*u_new; 
            
        if(~mod(iter,5)) 
            disp(['Iteration: ', num2str(iter), ' normResidual: ',num2str(normResidual)]);
        end
        iter=iter+1;
    end
    if(maxiter==iter)
        error('Failed to converge');
    end

    time1=toc(time1);
    disp(['Time to solve: ',num2str(time1),'s With ', num2str(iter),' iterations'])
    % Decomposition of the different terms
    sol(DOFl_v) = d(1:size(Kll,1));
    sol(DOFr_v) = dR_v;
    
    pl = d(size(Kll,1)+1:end);
    p=zeros(size(COOR_p,1),1);
    p(DOFl_p)=pl;
    p(DOFr_p)=dR_p;
    u = sol(1:2:end-1);
    v = sol(2:2:end);
end