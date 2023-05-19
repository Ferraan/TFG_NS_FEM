function  [Ne, BeXi] = Linear3N(xi)  
% Shape functions and derivatives for 3-node linear element 
 
% Matrix of shape functions
Ne =[-0.5*xi*(1-xi)  1-xi^2  0.5*xi*(1+xi)]; 
% Matrix of the gradient of shape functions 
BeXi = [-0.5+xi -2*xi 0.5+xi] ; 
end