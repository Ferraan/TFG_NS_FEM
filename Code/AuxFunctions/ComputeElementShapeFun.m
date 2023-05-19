function [weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% This function returns, for each "TypeElement" and 'TypeIntegrand'
% (K/RHS)*
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)
%
% *) TypeIntegrand = 'K': For integration of Conductance Matrix
% *) TypeIntegrand = 'RHS': For integration of flux vectors

Q2=1;
switch TypeElement
    case 'Linear'
        if nnodeE ==2
            [weig,posgp,shapef,dershapef] = Linear2NInPoints(TypeIntegrand) ;
        elseif nnodeE==3
            [weig,posgp,shapef,dershapef] = Linear3NInPoints(TypeIntegrand) ;
        else
            error('Option not implemented')
        end
    case 'Triangle'
        if nnodeE ==3
            [weig,posgp,shapef,dershapef] = Triangle3NInPoints ;
        else
            error('Option not implemented')
        end
    case 'Quadrilateral'
        if nnodeE ==4
            if (Q2)
              [weig,posgp,shapef,dershapef] = Quadrilateral9NInPointsPressure ;     
            else
               [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints ;
            end
        elseif nnodeE == 9
            [weig,posgp,shapef,dershapef] = Quadrilateral9NInPoints ;
        end
    otherwise
        error('Option not implemented')
end