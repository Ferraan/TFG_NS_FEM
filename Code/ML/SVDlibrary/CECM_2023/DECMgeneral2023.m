function [ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral2023(A,U,wFE,xFE,DATA_ECM,MESH,DATA)

if nargin  == 0
    load('tmp1.mat')
end

% List of elements to be excluded from the initial SET
DATA_ECM = DefaultField(DATA_ECM,'ListElementsToExclude',[]) ;

INDSEL = 1:length(wFE) ;

DATA_ECM = DefaultField(DATA_ECM,'IND_POINTS_CANDIDATES',INDSEL) ; 
 

if ~isempty(DATA_ECM.ListElementsToExclude)
    ListGaussToExclude = small2large(DATA_ECM.ListElementsToExclude,MESH.ngausE) ;
    INDSEL  = setdiff(DATA_ECM.IND_POINTS_CANDIDATES,ListGaussToExclude) ;
else
    INDSEL = DATA_ECM.IND_POINTS_CANDIDATES ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Empirical cubature method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_ECM = DefaultField(DATA_ECM,'errorDECM',0) ;
DATA_DECM = [] ;
DATA_DECM.TOL = DATA_ECM.errorDECM ;
DATA_DECM.IND_POINTS_CANDIDATES = INDSEL ;
DATA_DECM.STORE_INFO_ITERATIONS = 1; 
[zDECM,wDECM,~,DATAOUTdecm]= DiscreteEmpiricalCubatureMethod(U',wFE,DATA_DECM)  ;
ECMdata.xDECM = xFE(zDECM,:);
ECMdata.wDECM = wDECM ;
% Determining the indices of the associated elements
setElements = large2smallREP(zDECM,MESH.ngausE) ;
setElementsREP = large2smallINCLUDEREP(zDECM,MESH.ngausE) ;

disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
%clipboard('copy',num2str(setElements'));

HYPERREDUCED_VARIABLES.setPoints = zDECM ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElementsREP ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = wDECM ;  % Set associated WEights
HYPERREDUCED_VARIABLES.PHI = bsxfun(@times,U,1./sqrt(wFE)) ;


% ------------------------------------------------
% Computing integration errors due to the DECM
% -----------------------------------------------

DATA = DefaultField(DATA,'CalculateErrorIntegralSnapshotMatrices',1) ; % = 0 ;

if DATA.CalculateErrorIntegralSnapshotMatrices == 1
        ErrorApproximateIntegral2(A,HYPERREDUCED_VARIABLES.PHI,wFE,DATA_ECM,zDECM,wDECM,DATA)  ;
end
