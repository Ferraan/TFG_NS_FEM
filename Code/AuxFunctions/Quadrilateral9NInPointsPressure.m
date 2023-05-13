function [weig,posgp,shapef,dershapef] = Quadrilateral9NInPointsPressure
    % Computation of the information regarding the Gauss points, the shape
    % functions and the derivative of those shapefunctions for
    % 9Integrations points in an element of 4 points. 
    
    w1=5/9; w2=8/9; w3=5/9;
    weig= [w1*w1 w2*w1 w3*w1 w1*w2 w2*w2 w3*w2 w1*w3 w2*w3 w3*w3];
    p = sqrt(3/5);
    posgp=[-p   -p;
            0   -p;
            p   -p;
           -p    0;
            0    0;
            p    0;
           -p    p;
            0    p;
            p    p]';
        
    N = @(x,y) 0.25*[(1-x).*(1-y) (1+x).*(1-y) (1+x).*(1+y) (1-x).*(1+y)];
    B = @(x,y) 0.25*[-(1-y) (1-y) (1+y) -(1+y); -(1-x) -(1+x) (1+x) (1-x)];
    
    
    shapef = zeros(length(weig),4);
    for i=1:length(weig)
        shapef(i,:) = N(posgp(1,i),posgp(2,i));
    end
    
    dershapef = zeros(2,4,length(weig));
    
    for j=1:length(weig)
        dershapef(:,:,j) = B(posgp(1,j),posgp(2,j));
    end
end