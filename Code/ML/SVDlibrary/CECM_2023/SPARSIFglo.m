function [xGAUSS,wGAUSS,DATA_ECM,VAR_SMOOTH_FE,POLYINFO,AUXVAR] ...
    = SPARSIFglo(w,DATA_ECM,coorECM,VAR_SMOOTH_FE)
%dbstop('6')
if nargin ==0
    load('tmp1.mat')
    %  DATA_ECM.FORCE_POINTS_TO_REMAIN_INSIDE = 1 ;
    %  DATA_ECM.SECOND_STAGE_ITERATIONS = 10 ;
    %  DATA_ECM.ONE_STAGE_ITERATIONS = 0 ;
end

b=   DATA_ECM.ExactIntegral; % PHI'*W ; % Exact integral basis functions 
VAR_SMOOTH_FE.ExactIntegral = b;
% Starting Cubature rule
% ----------------------
wOLD =w ;
xOLD = coorECM ;
DATA_ECM.mINI = size(xOLD,1) ;
HISTORY = [] ; HISTORY.POINTS_all ={} ; 
HISTORY.WEIGHTS_all ={} ; 
HISTORY.POINTS_all{1} = xOLD ;
HISTORY.WEIGHTS_all{1} = wOLD ;
HISTORY.ISALLPOSITIVE = 1  ;


DATA_ECM = DefaultField(DATA_ECM,'THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS',5); %  Negative points allowed at each iteration
DATA_ECM = DefaultField(DATA_ECM,'SECOND_STAGE_ITERATIONS',0); % Refine iterations, 2nd stage
DATA_ECM = DefaultField(DATA_ECM,'maxITER_allowed_residual_withoutDECREASE',10); %  Number of iterations allowed without the residual experiencing any decrease

% List of design variables (DOFs ) 
% --------------------------------
[npoints, ~ ]= size(xOLD) ;
VARC.POINTSRp =[] ; % These are the indexes of the points for which   position   are constrained, but weights aren't. These conditions might
% vary during the iterations.  
VARC.POINTSl = 1:npoints; % These are the indexes of the points for which both weights and positions are unknowns
VARC.POINTSl(VARC.POINTSRp)  = [] ;
VARC.POINTSRpw = []; % These are the indexes of the points for which   weights and positions are constrained

% ---------------------------------
% FIRST STAGE ---direct elimination of points (by making its weights equal to zero) 
% ---------------------------------
DATA_ECM = DefaultField(DATA_ECM,'NumberOfCECMpoints',0) ; % The number of CECM points is given (MOD 8-April-2023) 

[xFIRSTS,wFIRSTS,DATA_ECM,POLYINFO,VARC,HISTORY] = SPARSIF_1stStage(xOLD,wOLD,b,DATA_ECM,VAR_SMOOTH_FE,VARC,HISTORY) ;

% --------------------------------------------------
% SECOND STAGE ---progressive elimination of points
% --------------------------------------------------
if DATA_ECM.SECOND_STAGE_ITERATIONS >1 && length(find(wFIRSTS>0)) >1 % && length(find(wFIRSTS>0)) >DATA_ECM.NumberOfCECMpoints
    [xSECONDS,wSECONDS,DATA_ECM,POLYINFO,VARC,HISTORY] = SPARSIF(xFIRSTS,wFIRSTS,b,DATA_ECM,VAR_SMOOTH_FE,VARC,HISTORY,POLYINFO) ;
end
% end

% Iterations with positive weights
% ----------------------------------
IterPositive= find(HISTORY.ISALLPOSITIVE ==1) ;
LastIteration = IterPositive(end) ;
% WE take this iteration as the final set of points
xNEWall = HISTORY.POINTS{LastIteration} ;
wNEWall = HISTORY.WEIGHTS{LastIteration} ;
INNZ = find(wNEWall>0) ; 
wNEW = wNEWall(INNZ) ; 
xNEW = xNEWall(INNZ,:) ; 
POLYINFO.setElements = POLYINFO.setElements(INNZ) ; 
% HISTORY.POINTS =  HISTORY.POINTS(1:LastIteration)  ;
% HISTORY.WEIGHTS =  HISTORY.WEIGHTS(1:LastIteration)  ;
HISTORY = DefaultField(HISTORY,'INDEXES_FIRST_STAGE',[]) ; 
HISTORY = DefaultField(HISTORY,'INDEXEX_SECOND_STAGE',[]) ; 

if  isempty(HISTORY.INDEXES_FIRST_STAGE)
    DATA_ECM.Include2ndStageIterations_PlotEvolutionWeights = 0  ; 
end
if DATA_ECM.Include2ndStageIterations_PlotEvolutionWeights ==1 
    IndSecondStageTake = LastIteration-length(HISTORY.INDEXES_FIRST_STAGE) ;
    if IndSecondStageTake <=0
        HISTORY.INDEXES_FIRST_STAGE = HISTORY.INDEXES_FIRST_STAGE(1:LastIteration) ;
        HISTORY.INDEXEX_SECOND_STAGE = {} ;
        HISTORY.POINTS_all =  HISTORY.POINTS_all(1:LastIteration)  ;
    HISTORY.WEIGHTS_all =  HISTORY.WEIGHTS_all(1:LastIteration)  ;
    else
        HISTORY.INDEXEX_SECOND_STAGE = HISTORY.INDEXEX_SECOND_STAGE(1:IndSecondStageTake) ;
        LastIndex = HISTORY.INDEXEX_SECOND_STAGE{end}(end) ; 
        HISTORY.POINTS_all =  HISTORY.POINTS_all(1:LastIndex)  ;
    HISTORY.WEIGHTS_all =  HISTORY.WEIGHTS_all(1:LastIndex)  ;        
    end
else     
    HISTORY.POINTS_all =  HISTORY.POINTS_all(1:LastIteration)  ;
    HISTORY.WEIGHTS_all =  HISTORY.WEIGHTS_all(1:LastIteration)  ;
end

%HISTORY.ELEMENTS_CONTAINING_POINTS =  HISTORY.ELEMENTS_CONTAINING_POINTS(1:LastIteration)  ;

disp('---------------------------------')
disp('Summary')
disp('-------------------------------')
disp(['Reduction in first stage: from ',num2str(length(wFIRSTS)),' to ',num2str(length(find(wFIRSTS ~=0)))])

disp('--------------------------------------------------')
disp(['Final integration rule with m =',num2str(length(wNEW)),' POINTS  (of ',num2str(size(coorECM,1)),'). ',...
    '. Rank Basis = ',num2str(length(b))]);
disp('--------------------------------------------------') ;
disp('Integration error') ;
disp('-------------------------------------------------------')
[PHIk_y,~,POLYINFO]=     EVALBASIS(xNEW,DATA_ECM,VAR_SMOOTH_FE,POLYINFO) ;
bNEW = PHIk_y'*wNEW ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error (%) =',num2str(errorINT)]) ;
%---------------
wGAUSS = wNEW ;
xGAUSS = xNEW ;
%
AUXVAR.HISTORY = HISTORY; 


 

