function [ECMdata,DATA_OUT]= ContinuousECMgen2023(A,xFE,wFE,DATA,MESH,DATA_ECM)
% Adaptation of ContinuousECM.m, general functions
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/08_BACK_TO_POLY/README_08.mlx
if nargin == 0
    load('tmp1.mat')
end
ECMdata = [] ;
%DATA_ECM.NameFileMesh_ECM = [DATA.BASE_FOLDER,'DECMpoints'];
%DATA_ECM.NameFileMesh_CECM = [DATA.BASE_FOLDER,'CECMpoints'] ;
% Basis matrix for the column space of A
% ------------------------------------------------------------------
DATA_ECM = DefaultField(DATA_ECM,'Method_Evaluation_Basis_Integrand',1) ;  % % 1 --> FE fitting   / 2 ---> FE analytical evaluation
DATA_ECM = DefaultField(DATA_ECM,'NumberOfCECMpoints',[]) ;  % % 1 --> 
if    isempty(DATA_ECM.NumberOfCECMpoints)
    DATA_ECM.NumberOfCECMpoints = 0  ;  
end

DATA = DefaultField(DATA,'ExactIntegral',[]) ;
if isempty(DATA.ExactIntegral)
    error('The exact integral should be calculated a priori')
end
[BasisA_squareW,VSinv] = GetBasisMatrixECM2(A,wFE,DATA_ECM) ;
%% DISCRETE EMPIRICAL CUBATURE METHOD (DECM)
% -------------------------------------------------------------------------
%[ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral(A,BasisA_squareW,wFE,xFE,DATA_ECM,MESH,DATA) ;
[ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral2023(A,BasisA_squareW,wFE,xFE,DATA_ECM,MESH,DATA) ;

% PLOT ORTHOGONAL FUNCTIONS TO BE INTEGRATED
DATA_ECM = DefaultField(DATA_ECM,'PLOT_ORTHOGONAL_FUNCTIONS_TO_BE_INTEGRATED',1) ;
PlotOrthogonalFunctToIntegrateCECM(HYPERREDUCED_VARIABLES.PHI,MESH,DATA_ECM) ; 


ECMdata.HISTORY.xDECM = DATAOUTdecm.HistoryPoints ;
ECMdata.HISTORY.wDECM = DATAOUTdecm.HistoryWeights ;
% ---------------------------------------------------------------------------
% DATA PREPARATION FOR Continuous EMPIRICAL CUBATURE METHOD (CECM)
% ---------------------------
DATA_ECM = DefaultField(DATA_ECM,'UseDECMpoints',0) ;
 
if DATA_ECM.UseDECMpoints == 0
    [VAR_SMOOTH_FE,DATA,DATA_ECM] = DataPreparationCECM(MESH,DATA,DATA_ECM,wFE,xFE,VSinv) ;
    
    % CONTINUOS EMPIRICAL CUBATURE METHOD (CECM)
    % -----------------------------------------------------------------------------
    [xNEW,wNEW,POLYINFO ,VAR_SMOOTH_FE,Ninterpolation,DATA_ECM,AUXVAR ] = ContinuousECM2023...
        (wFE,xFE,HYPERREDUCED_VARIABLES,DATA_ECM,VAR_SMOOTH_FE)  ;
    % --------------------------------------------------------------------------------------------
    
    %%%  PREPARING OUTPUT/PLOTTING FACILITIES
    % --------------------------------------------
    ECMdata.xCECM = xNEW ;
    ECMdata.wCECM = wNEW ;
    ECMdata.setElements = POLYINFO.setElements ; %  ELEMENTS_xNEW ;
    ECMdata.Ninterpolation = Ninterpolation ;
    ECMdata.INTexac = DATA.ExactIntegral;
    AUXVAR.DATAOUTdecm = DATAOUTdecm;
    
    VariousPLOTS_CECM2023(DATA_ECM,ECMdata,MESH,AUXVAR,VAR_SMOOTH_FE)  ;
    
    DATA_OUT.DATA_ECM = DATA_ECM ;
    DATA_OUT.ECMdata = ECMdata ;
    DATA_OUT.MESH = MESH ;
    DATA_OUT.AUXVAR = AUXVAR ;
    DATA_OUT.VAR_SMOOTH_FE = VAR_SMOOTH_FE ;
    DATA_OUT.DATAOUTdecm = DATAOUTdecm ;
    DATA_OUT.POLYINFO = POLYINFO ;
    
    DATA_OUT.indexPoints_DECM = HYPERREDUCED_VARIABLES.setPoints;
    
else
    
    ECMdata.xCECM =[];
    ECMdata.wCECM = [];
    ECMdata.setElements = HYPERREDUCED_VARIABLES.setElements ; %  ELEMENTS_xNEW ;
    ECMdata.Ninterpolation = [] ;
  %  ECMdata.INTexac = DATA.ExactIntegral;
  %  AUXVAR.DATAOUTdecm = DATAOUTdecm;
    
  %  VariousPLOTS_CECM2023(DATA_ECM,ECMdata,MESH,AUXVAR,VAR_SMOOTH_FE)  ;
    
    DATA_OUT.DATA_ECM = DATA_ECM ;
    DATA_OUT.ECMdata = ECMdata ;
    DATA_OUT.MESH = MESH ;
    DATA_OUT.AUXVAR = [] ;
    DATA_OUT.VAR_SMOOTH_FE = [] ;
    DATA_OUT.DATAOUTdecm = [] ;
    DATA_OUT.POLYINFO = [] ;
    
    DATA_OUT.indexPoints_DECM = HYPERREDUCED_VARIABLES.setPoints;
    
end




