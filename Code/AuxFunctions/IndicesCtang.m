function   [i,j] = IndicesCtang(m,p)

if nargin == 0
    p = 2;
    m = 6 ;
end

n = m ;
ngaus = m/p ;
nzmax = m*p ;
ij = zeros(p^2,2) ;
for ielem = 1:p
    for jelem = 1:p
        iniELEM = (ielem-1)*p +jelem;
        ij(iniELEM,1) = ielem ;
        ij(iniELEM,2) = jelem ;
    end
end

ij = repmat(ij,ngaus,1);
indSUM = 0:p:(ngaus-1)*p ;
indSUM = repmat(indSUM,p^2,1);
indSUM = reshape(indSUM,size(indSUM,1)*size(indSUM,2),1) ;
ij = bsxfun(@plus,ij,indSUM) ;
i = ij(:,1);
j = ij(:,2);