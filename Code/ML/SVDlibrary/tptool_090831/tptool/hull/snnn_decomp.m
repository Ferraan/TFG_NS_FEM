function [W V] = snnn_decomp(U)
%SNNN_DECOMP snnn decomposition of a matrix with orthogonal columns.
%	[W V] = SNNN_DECOMP(U)

% TODO: doc

sumtol = 1e-5;
Usum = sum(U, 1);
X = diag(Usum);
% check if any element of the row sums equals to zero
for i = find(abs(Usum) < sumtol)
	X(i,i) = X(i,i) + 1;
	if i == 1
		X(i,i+1) = -1;
	else
		X(i,i-1) = -1;
	end
end

if all(abs(sum(U*U',2)-1) < sumtol)
	% check if X is invertible
	if all(abs(Usum) < sumtol)
		disp('SN transformation is close to singular, expect errors');
	end
	Usn = U*X;
	%Vsn = inv(X);
else
	Usn = [U 1-sum(U,2)];
	% TODO: calc V on the fly?
	%ns = size(U, 2);
	%Vsn = eye(ns, ns + 1);
end

% >= 0 condition:
%   M = d*ones(n,n) + c*eye(n,n)
%   SN --> c + n*d == 1 --> det(M) == c^n + n*d*c^(n-1) == c^(n-1)
%   NN --> 1/(1 - umax n) <= c <= 1/(1 - umin n)
n = size(Usn, 2);
umin = min(min(min(Usn)), 0);
umax = max(max(max(Usn)), 2/n);
c1 = 1/(1 - umin*n);
c2 = 1/(1 - umax*n);
if c1 < -c2
	c = c2;
else
	c = c1;
end
d = (1 - c)/n;
M = d * ones(n,n) + c * eye(n,n);

W = Usn*M;
if nargout==2
	V = W\U;
end
