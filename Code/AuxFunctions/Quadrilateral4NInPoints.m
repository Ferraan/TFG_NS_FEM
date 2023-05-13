function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints 
    % Computation of the information regarding the Gauss points, the shape
    % functions and the derivative of those shapefunctions for
    % 4Integrations points in an element of 4 points. 
    
    weig = [1 1 1 1];
    posgp = 1/sqrt(3)*[-1 1 1 -1; 
                       -1 -1 1 1];
    
    N = @(x,y) 0.25*[(1-x).*(1-y) (1+x).*(1-y) (1+x).*(1+y) (1-x).*(1+y)];
    B = @(x,y) 0.25*[-(1-y) (1-y) (1+y) -(1+y); -(1-x) -(1+x) (1+x) (1-x)];
    shapef = zeros(length(weig),length(weig));
    for i=1:length(weig)
        shapef(i,:) = N(posgp(1,i),posgp(2,i));
    end
    
    dershapef = zeros(2,length(weig),length(weig));
    
    for j=1:length(weig)
        dershapef(:,:,j) = B(posgp(1,j),posgp(2,j));
    end
end




