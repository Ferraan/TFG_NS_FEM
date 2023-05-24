%az U1 ponthalmazt burkolja polarv alaku tetra�derrel
%fi2 a tetra�der cs�csai
% U1=U2*fi2

function [U2,v]=polarorto(U1,polarv)

[n r]=size(U1);

a=polarv;

for i=1:r-2    
    v(:,i)=sin(a(:,i)).*prod(cos(a(:,1:i-1)),2);
end

v(:,r-1)=prod(cos(a(:,:)),2);
v(:,r)=0;
v(:,r)=1-sum(v')';

v=v+ones(r,1)*(U1(1,:)-sum(v)/r);  %eltolja a szimplexet, hogy az U1(1,:) pont a belsej�ben legyen

for i=1:r
    U2=U1*inv(v);
    m=min(U2(:,i));
    v=v-m*(v-ones(r,1)*v(i,:)); %i. pontb�l nagy�t, h. a szemk�zti lap �rintse a ponthslmazt
    v(:,r)=0;
    v(:,r)=1-sum(v')';   
end

U2=U1*inv(v);
v;