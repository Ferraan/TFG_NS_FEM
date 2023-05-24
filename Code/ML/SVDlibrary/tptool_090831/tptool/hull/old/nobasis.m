function [U, V] = nobasis(A, ia)
%NOBASIS Column space basis with SNNN and NO properties.
%	[U, V] = NOBASIS(A, inter_active)

[U, V] = nobasis(A, ia);
[sor, dim] = size(U);

if dim==2
	P = [1 0; -1 1];
	XM = U/P;
	mine = min(XM(:,2));
	maxe = max(XM(:,2));
	TG = [1 mine; 1 maxe];
elseif dim==3
	P = [1 0 0; -1 1 0; 0 -1 1];
	XM = U/P;
	x = XM(:,2); %x=XM(1:10,2);
	y = XM(:,3); %y=XM(1:10,3);
	figure(2);
	plot(x,y,'x');
	hold on;
	G = XM;
	n1 = size(U, 1);

	Igo = input('input: 1->input 3 row #, 2->Geometric Geometry , 3->input 3 points, ');
	if Igo==1
		Imin = input('Input convex Hull #1:  ');
		Imed = input('Input convex Hull #2:  ');
		Imax = input('Input convex Hull #3:  ');
		TG = G([Imin Imed Imax],:);
	elseif Igo==2
		m1 = (G(2,3)-G(1,3))/(G(2,2)-G(1,2));
		m2 = (G(n,3)-G(n-1,3))/(G(n,2)-G(n-1,2));
		x = (m1*G(1,2)-m2*G(n-1,2)-(G(1,3)-G(n-1,3)))/(m1-m2);
		y = m1*(x-G(1,2))+G(1,3);
		TG = [G(1,:);G(n,:);1 x y];
	elseif Igo==3
		Imin = input('Input point #1:  ');
		Imed = input('Input point #2:  ');
		Imax = input('Input point #3:  ');
		TG = [1 Imin;1 Imed;1 Imax];
	end
	figure(2)
	TTT = [TG;TG(1,:)];
	% subplot(3,2,1)
	plot(TTT(:,2),TTT(:,3),'g')     
	hold off;
else
	error('dimension error: dim = %d', dim);
end
U = XM/TG;
V = TG*P*V;
%h=max(max(abs(Y-U*V)))
