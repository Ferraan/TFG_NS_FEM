function [DOFl_v, DOFr_v, uBAR] = BC_LidDriven(COOR_v,CNb_v,rnod_v,ndim)
% This function is meant to compute the restricted and free degrees of freedom
% as well as its corresponding vector of restricted velocities.
% Inputs:
%   - COOR_v: Matrix containing the coordinates of the nodes for the v-mesh
%   - CNb_v: Matrix containing the connectivities at the boundaries for
%   v-mesh.
%   - rnod_v: Vector containing all the restricted nodes of the velocities
%   - ndim: Should be 2 (3D not implemented yet)
% Outputs:
%   - DOFl_v: vector with the free desgrees of freedom for v
%   - DOFr-v: vector with the restricted degrees of freedom for v
%   - uBAR: Vector containing the preseribed velocities for the restricted
%   degrees of freedom for v

global Q2;

    % Free and restricted degrees of freedom
    DOFl_v = (1:ndim*size(COOR_v,1))';
    DOFr_v = Nod2DOF(rnod_v,ndim);
    DOFl_v(DOFr_v) = [];

    % Computation of restricted degrees of freedom for the velocity case
    uBAR = zeros(size(COOR_v,1)*ndim,1);
    for i = 1:size(CNb_v,1)
        if(COOR_v(CNb_v(i,1),2) == 1)
            uBAR(Nod2DOF(CNb_v(i,1),ndim)) = [1 0]'; %Case of beeing at the lid
        end
        if (Q2) 
            if (COOR_v(CNb_v(i,3),2) == 1)
                uBAR(Nod2DOF(CNb_v(i,3),ndim)) = [1 0]';%In case that the mesh is Q2, this other node shall also be checked
            end
        end
    end
    uBAR(DOFl_v) = [];
end