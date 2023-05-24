function  [xOLD,wOLD,DATALOC,POLYINFO,VARC,HISTORY] =...
    SPARSIF(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,POLYINFO)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp1.mat')
end

CONVERGENCE = 1 ;
iter = 1;
VARC.POINTSl = find(wOLD ~=0);% These are the indexes of the points for which both weights and positions are unknowns
% (candidates for being eliminated)
VARC.POINTSRp = []; % These are the indexes of the points for which   positions   are constrained, but weights aren't
VARC.POINTSRpw = find(wOLD ==0); % These are the indexes of the points for which   weights and positions are constrained
DATALOC = DefaultField(DATALOC,'Include2ndStageIterations_PlotEvolutionWeights',1) ; % Plotting option
DATALOC = DefaultField(DATALOC,'VARIABLE_ITERATIONS_SECOND_STAGE',0) ; % Number of iterations is variable, see Eliminate1PointPROG.m

HISTORY.INDEXES_FIRST_STAGE = 1:length(HISTORY.POINTS_all) ; % For Plotting purposes
HISTORY.INDEXEX_SECOND_STAGE = {};% For Plotting purposes
HISTORY.CONTROL_POINTS = [] ; 
HISTORY.POINTS = HISTORY.POINTS_all ; 
HISTORY.WEIGHTS = HISTORY.WEIGHTS_all ; 

%HISTORY.CONTROL_POINTS =[] ;
NumberIterFirstStage = length(HISTORY.INDEXES_FIRST_STAGE) ;
CurrentITERglo = NumberIterFirstStage  ;
% ------------------------------------------------------
while CONVERGENCE ==1  &&  length(find(wOLD>0)) >DATALOC.NumberOfCECMpoints
    DATALOC.iter = iter ;
    %  [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC] = ControlPointsAlgLARGE(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC) ;
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC ] = MAKE1ZERO(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,POLYINFO) ;
    if CONVERGENCE == 1
        
        [HISTORY,INDnneg,CurrentITERglo] = HistoryPointsUpdate2023(wNEW,xNEW,HISTORY,DATALOC,iter,CurrentITERglo) ;
        DATALOC.HISTORY_LOCAL.x = {} ;
        DATALOC.HISTORY_LOCAL.w  = {} ;
        iter = iter + 1 ;
        xOLD = xNEW ;        wOLD = wNEW ;
        
        if   length(find(wOLD>0))== DATALOC.NumberOfCECMpoints
            disp('Number of points equal to the number prescribed to the user. ..Exiting')
            break
        end
        
    end
end
%HISTORY.CONTROL_POINTS = DATALOC.HISTORY_LOCAL.ControlPoints ;

if iter == 1
    INDnneg = find(wOLD ~=0);
end

disp('------------------------------------------------------------------------------------------------')
disp(['SECOND STAGE: Integration rule with ',num2str(length(INDnneg)),' of ',num2str(length(wOLD)),' points'])
disp('------------------------------------------------------------------------------------------------')
 
