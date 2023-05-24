function [ip,cs] = ino(U)
% INO - m�trix felbont�s az INO felt�tel szerint
%
% [IP,BASE]=INO(U)
% 
% Az U m�trix nxr-es, n>r �s SN, azaz a sorok �sszege
% 1 kell legyen;
% 
% A kimenet az nxr-es IP �s az rxr-es BASE m�trix, ahol
% a BASE SN, az IP pedig SN, NN(nemnegat�v elemu) 
% �s INO(minden oszlop legkisebb eleme 0
% 
% �s amelyekre U=IP*BASE

%U=[.5 .1 .2 .2;.1 .1 .1 .7;.2 0 .4 .4; .2 .8 0 0;.3 .3 .2 .2];
[n r] = size(U);
UU = U;
U(:,r) = 0;

egy = min(U);  % a burkol�szimplex lapjainak egyenlete/1
egy(r) = max(sum(U,2));
cs = ones(r,1)*egy(1:r-1); % a burkol�szimplex cs�csai/1
for a = 1:r-1
    cs(a,a) = 2*egy(r) + egy(a) - sum(egy); % a burkol�szimplexcsucsai/2
end
cs(:,r) = 1 - sum(cs,2); %az r. dimenzi�
%for a=1:r
%    RNO(a,:)=sum(cs)-(r-1)*cs(a,:);% av�gleges szimplex cs�csai
%end

%ip=UU*inv(RNO); %az interpol�ci�s m�trix
ip = UU*inv(cs);  %az interpol�ci�s m�trix
