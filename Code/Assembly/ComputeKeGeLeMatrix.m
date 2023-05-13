function [Ke, Ge, Le] = ComputeKeGeLeMatrix(weig,~,dershapef_v,shapef_p,dershapef_p,Xe_v,Xe_p,nu)
% This function is meant to compute the element K, G and L matrices
% Inputs:
%   - weig: Vector with the weights of the Gauss integration
%   - shapef_v: Shape functions of velocity evaluated at the Gauss points
%   - dershapef_v: Derivative of shape functions (vel) evaluated at Gauss points
%   - shapef_p: Shape functions of velocity evaluated at the Gauss points
%   - dershapef_p: Derivative of shape functions (p) evaluated at Gauss points
%   - Xe_v: matrix with the coordinates of the nodes (vel) in the original
%   domain
%   - Xe_p: matrix with the coordiantes of the nodes (p) in the original domain
% Outputs:
%   - Ke: Element conductance matrix
%   - Ge: Element gradient matrix
%   - Le: Element GLS stabilisation matrix


nu;
global a0;


% Parameters for the problem
ndim = size(Xe_v,1) ; ngaus = length(weig) ; nnodeE_v = size(Xe_v,2)  ;
nnodeE_p = size(Xe_p,2);

% Initialization of element matrices
Ke = zeros(nnodeE_v*ndim,nnodeE_v*ndim);
Ge = zeros(nnodeE_p,nnodeE_v*ndim);
Le = zeros(nnodeE_p,nnodeE_p);

% Parameters for the GLS stabilisation
v1 = Xe_p(:,3)-Xe_p(:,2);
v2 = Xe_p(:,1)-Xe_p(:,2);
he = norm(cross([v1;0],[v2;0]));

tau1 = a0*(he^2/4*nu);

for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi_v = dershapef_v(:,:,g) ; 
    BeXi_p = dershapef_p(:,:,g) ; 

    % Jacobian Matrix 
    Je_v = Xe_v*BeXi_v' ;
    Je_p = Xe_p*BeXi_p' ; 
    
    % Jacobian 
    detJe = det(Je_v) ; 
    % Matrix of derivatives with respect to physical coordinates 
    BeTILDE_v = inv(Je_v)'*BeXi_v ;
    BeTILDE_p = inv(Je_p)'*BeXi_p ;
    
    % Matrix of symmetric gradient 
    Be_v = QtransfB_special(BeTILDE_v,ndim);
    % Element matrices
    Ke = Ke + weig(g)*detJe*(Be_v'*nu*Be_v);
    Ge = Ge +weig(g)*detJe*(shapef_p(g,:)'*[1 0 0 1]*Be_v);% 
    if (a0 ~=0 )
        Le = Le + weig(g)*detJe*tau1*(BeTILDE_p'*BeTILDE_p);
    end
end