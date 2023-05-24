function  [xGOOD,wGOOD,DATALOC,POLYINFO,VARC,HISTORY] = SPARSIF_1stStage(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY)
% First stage of reduction (direct removal of weights)
% .....................................................
if nargin == 0
    load('tmp.mat')
end

CONVERGENCE = 1 ;
iter = 1;
wNEW = [] ;
INDnneg = 1:length(wOLD) ;
xGOOD = xOLD ; wGOOD = wOLD ;
POLYINFO.setElements = VAR_SMOOTH_FE.setElements; %  Set of elements containing xNEW (initial guess)
POLYINFO.TriangulationDelaunay = cell(size(VAR_SMOOTH_FE.CN,1),1) ;   % Triang. For Delaunay 3D search 
DATALOC = DefaultField(DATALOC,'TOL_NewtonRaphson_EliminationPoints',1e-8) ;
DATALOC = DefaultField(DATALOC,'MaxIterationsNR_ElimPoints',40) ;

% ------------------------------------------------------
while CONVERGENCE ==1  %&& length(xOLD)>1
    DATALOC.iter = iter ;
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC] = MAKE1ZERO_1STEP(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,POLYINFO) ;
    if CONVERGENCE == 1
        INDnneg = find(wNEW ~=0);
        HISTORY.POINTS_all{end+1} = xNEW ;
        HISTORY.WEIGHTS_all{end+1} = wNEW  ;        
        ALLPOSI = all(wNEW(INDnneg)>0)  ;
        if ALLPOSI
            xGOOD = xNEW ;
            wGOOD = wNEW ;
        end
        HISTORY.ISALLPOSITIVE(end+1) = ALLPOSI ;
        iter = iter + 1 ;
        xOLD = xNEW ;        wOLD = wNEW ;
        if length(INDnneg) == 1 
            break
        end
        
        if length(INDnneg) ==   DATALOC.NumberOfCECMpoints
            disp(['Number of CECM points equal to the number prescribed by the user...'])
            break
        end
        
    end
end

if ~isempty(wNEW)
    disp('------------------------------------------------------------------------------------------------')
    disp(['FIRST STAGE: Integration rule with ',num2str(length(INDnneg)),' of ',num2str(length(wNEW)),' points'])
    disp('------------------------------------------------------------------------------------------------')
else
    error(['Convergence not achieved; try to decrease the tolerancd for the Newton-Raphson'])
end

