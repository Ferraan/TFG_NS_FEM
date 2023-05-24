function [U,S,V,eSVD] = SVDT(B,epsilon,DATA)
%  SVDT(B,epsilon) computes a truncated SVD of B (with   tolerance epsilon)
%  It employs the built-in matlab function "svd", but it only returns the
% singular values  truncated according to the   the specified tolerance
%
% S = SVDT(B,epsilon) gives the truncated vector of singular values
%
% [U,S] = SVDT(B,epsilon) furnishes the truncated matrix of left singular vectors and
% the associated singular values
%
% [U,S,V] = SVDT(B,epsilon)  returns the full factorization (truncated).
%
% By default,  norm(U*diag(S)*V'-B,'fro') <= epsilon
% If epsilon = 0, then the truncation criterion is the epsilon machine
% parameter
% If the third argument REL = 1, (by default REL = 0), then
% norm(U*diag(S)*V'-B,'fro') <= epsilon*norm(B,'fro')*
%
% -----------------------------------
% Written by Joaquín  A. Hernández, March 2017.
% UPC/CIMNE, Barcelona, Spain
% jhortega@cimne.upc.edu
if nargin == 0
    load('tmp1.mat')
    epsilon = 0 ;  DATA = [] ;
end
if nargin ==1
    epsilon = 0 ;  DATA = [] ;
elseif nargin ==2
    DATA = [] ;
end

DATA = DefaultField(DATA,'COMPUTE_U',1) ;
DATA = DefaultField(DATA,'COMPUTE_V',1) ;
DATA = DefaultField(DATA,'RELATIVE_SVD',0) ;
DATA = DefaultField(DATA,'SVD',[]) ;
DATA.SVD = DefaultField(DATA.SVD,'MaxSize',max(size(B))) ;



CALCU =DATA.COMPUTE_U ; CALCV = DATA.COMPUTE_V ; U = [] ; V= [] ;
dimMATRIX = DATA.SVD.MaxSize;
M = size(B,1); N = size(B,2);
if M>=N
    if CALCU==1 & CALCV==1
        [U, S,V] = svd(B,'econ') ;  % U --> M xN, V --> N x N
    elseif CALCU==1 & CALCV==0
        [U, S] = svd(B,'econ') ;  % U --> M xN, V --> N x N
    else
        S = svd(B,'econ') ;
    end
else
    % If N>M, it proves more efficient to perform the SVD of B'
    if CALCU==1 & CALCV==1
        [V, S,U] = svd(B','econ') ; % U --> M x M, V --> Nx M
    elseif  CALCU==1 & CALCV==0
        [~, S,U] = svd(B','econ') ; % U --> M x M, V --> Nx M
    else
        S = svd(B','econ') ;
    end
end
if size(S,2)>1
    S = diag(S) ; % S = [S1 S2 ... SR]'
end

tol = dimMATRIX*eps(max(S));

R = sum(S > tol);  % Definition of numerical rank


eSVD = 0 ;

%dbstop('73')
if  epsilon ==0
    K = R ;
else
    %    dbstop('70')
    SingVsq =  (S.*S) ;
    SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
    normEf2 = sqrt(cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
    
    if DATA.RELATIVE_SVD ==1
        epsilon = epsilon*normEf2(end) ;
        
    end
    T = (sum(normEf2<epsilon)) ;
    K = length(S)-T ;
    
    
end
K = min(R,K);

%dbstop('88')
if length(S) >K
    eSVD = sqrt(sum(S(K+1:end).^2)) ;
end

S = S(1:K);
if  ~isempty(U) & isempty(V)
    U = U(:,1:K)  ;
    %  varargout={U,S};
elseif ~isempty(U) & ~isempty(V)
    U = U(:,1:K)  ;
    V =  V(:,1:K) ;
    % %
    % else
    %     varargout = {S} ;
end




