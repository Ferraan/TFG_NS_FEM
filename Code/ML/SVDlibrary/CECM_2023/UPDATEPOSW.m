function [xNEW, wNEW, nF,POLYINFO,ISOUT,VARCnew  ] = UPDATEPOSW(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp.mat')
end
% Evaluation of the function and the gradients at all points%%
% --------------------------------------------------------
%setElementsBefore = POLYINFO.setElements; 
[PHIk_y,dPHIk_y,POLYINFO]=     EVALBASIS(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO);
% ------------------------------
bNEW = PHIk_y'*wNEW ;
ndim = size(xNEW,2) ;
% Residual
r = b-bNEW ;
nF = norm(r)  ;
% Computation of the Jacobian matrix (derivative with respect to X)
J = zeros(length(r),length(wNEW)*ndim);
for idim = 1:length(dPHIk_y)
    J(:,idim:ndim:end) =  bsxfun(@times,dPHIk_y{idim},wNEW)' ;
end
% We are only interested in the components associated to POINTSl
ndim = size(xNEW,2) ;
DOFl = small2large(VARCnew.POINTSl,ndim) ;
VARCnew.DOFl  = DOFl ;
J = J(:,DOFl) ;
% The other two blocks are formed by
POINTS_F = [VARCnew.POINTSl(:);VARCnew.POINTSRp(:) ] ;
J = [J,  PHIk_y(POINTS_F,:)'] ;
DATALOCSVD.RELATIVE_SVD = 1;
TOL = 1e-10 ;
[UU,SS,VV] = SVDT(J,TOL,DATALOCSVD) ;
% Solution of the underdetermined system of equations (more unknowns than equation)
if length(SS) == size(J,1)
    % No need to correct rank
    % -------------------------------------------------------
    delta_q = J\r ;    % Matlab produces a sparse solution (minimize  the l1 norm of delta_q)
    POLYINFO.MESSAGE_RANK = '' ;
else
   % disp(['Incomplete rank: number of equations = ',num2str(size(J,1)),'; rank = ',num2str(length(SS)),' (TOLrel =',num2str(TOL)  ,')']);
    UU = bsxfun(@times,UU',1./SS)'  ;
    delta_q = VV'\(UU'*r) ;
    POLYINFO.MESSAGE_RANK = ['(Rank ',num2str(length(SS)),' of ',num2str(size(J,1)),')'] ; 
end

[xNEW,wNEW,ISOUT,VARCnew,POLYINFO] =  UpdateAndCheckOutside2023(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;

