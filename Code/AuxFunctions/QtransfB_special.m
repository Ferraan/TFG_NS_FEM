function Be = QtransfB_special(BeTILDE,ndim)
% This function is meant to distribute the derivatives of the shape
% function as required to create the gradient operator for fluid problems.

if nargin == 0
    load('tmp2.mat')
end


nnodeE = size(BeTILDE,2); 
if ndim==2 
   nstrain = 3 ;    Be = zeros(nstrain,nnodeE*ndim) ; 
   column1 = 1:2:(nnodeE*2-1) ;
   column2 = 2:2:nnodeE*2 ;
   Be(1,column1) = BeTILDE(1,:) ; 
   Be(2,column2) = BeTILDE(1,:) ; 
   Be(3,column1) = BeTILDE(2,:) ; 
   Be(4,column2) = BeTILDE(2,:) ; 
elseif ndim ==3 
    error('Not implemented, sorry')
else
   error('Incorrect option')
end