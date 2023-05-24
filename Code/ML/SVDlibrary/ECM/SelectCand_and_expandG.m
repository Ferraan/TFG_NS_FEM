function [G,b,sqW,DATA,y,yORIG,TOL,Gnorm,z] = SelectCand_and_expandG(G,W,DATA)

% BasisF is already multiplied by sqrt(W)
% -----------------------------------------
% BasisF = G'
% Orthogonal complement of sqrt(W)
% See
% Expansion with the orthogonal complement
% ----------------------------------------
sqW = sqrt(W) ;
a = sqW - G'*(G*sqW ) ;
V =sum(W) ;
aN = norm(a) ;
importanceA = aN^2/V ;
TOLLOC = 1e-10 ;
if importanceA > TOLLOC
    G = [ a'/aN;G] ;
end
% -------------------------------------------

b = G*sqW ;  % b Vector (exact integral)

Gnorm =sqrt(sum(G.*G,1)) ; % Norm of Columns of G
M = size(G,2) ;  % Number of FE points
DATA = DefaultField(DATA,'TOL',0) ; % Default tolerance for convergence
TOL = DATA.TOL ;
% INITIALIZATIONS
% ------------------------
z = [] ; % Set of integration points
y=1:M ;
DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-6) ;
% GnormNOONE =sqrt(sum(G(1:end-1,:).*G(1:end-1,:),1)) ; % Norm of Columns of G(see old, wrong line above)
%GnormNOONE =sqrt(sum(G(2:end,:).*G(2:end,:),1)) ; % cHANGE, 23-mARCH-2022 
if DATA.TOLFilterCandidatePoints >0
    TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
    %rmvpin = find(GnormNOONE(y)<TOL_REMOVE) ;
    rmvpin = find(Gnorm(y)<TOL_REMOVE) ;
    y(rmvpin) = [] ;
end
DATA = DefaultField(DATA,'RemoveColumnsWithNegativeProjection',0); % 4-Dec-2019

DATA = DefaultField(DATA,'IND_POINTS_CANDIDATES',[]) ;

if ~isempty(DATA.IND_POINTS_CANDIDATES)
    y = intersect(y,DATA.IND_POINTS_CANDIDATES) ;
end
yORIG = y ;

