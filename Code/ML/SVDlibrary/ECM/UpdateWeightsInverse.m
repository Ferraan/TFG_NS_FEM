function [Bast x] = UpdateWeightsInverse(A,Aast,a,xold,r)

if nargin == 0
    format long g
    P = 100 ; 
    k = 49 ; 
    A = randn(P,k) ; 
    a = randn(P,1) ;
    Aast = inv(A'*A) ; 
    B = [A a] ; 
    BastReal = inv(B'*B) ;
    b = randn(P,1) ; 
    xold = Aast*(A'*b) ; 
    r = b - A*xold ; 
    xREAL = BastReal*(B'*b) ;
    
 
end
c = A'*a ; 
d = Aast*c ; 
s = a'*a-c'*d ; 
Bast = [Aast + (d*d')/s  -d/s; -d'/s  1/s] ; 
v = a'*r/s ; 
x = [(xold -d*v ); v] ;

 
