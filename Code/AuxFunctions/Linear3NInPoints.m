function [weig,posgp,shapef,dershapef] = Linear3NInPoints(TypeIntegrand) ;
% This function returns, for each 2-node 1D linear element,
% and using a 1 points Gauss rule (ngaus=1): , if TypeIntegrand = K
% and ngaus=2, if TypeIntegrand = RHS:
%
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)

switch TypeIntegrand
    case {'K'}
        error('Option not implemented')
    case {'RHS'}
        % Two integration points
        weig  = [5/9 8/9 5/9] ;
        posgp = sqrt(0.6)*[-1 0 1 ] ;
end

ndim = 1; nnodeE = 3 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xi = posgp(g) ;
    [Ne BeXi] = Linear3N(xi) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end