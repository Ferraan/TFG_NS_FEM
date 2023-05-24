clear;

addpath("Postprocess\")

files=dir("GIDPOST/LidDriven*NSS*");
DATA=0;
i=1;
for file=files'
    load(strcat("GIDPOST\",file.name,'/',file.name,".mat"));
    U=SNAP_cluster.u.U(DOFl,:);
    GidPostProcessModes(COOR,CN,TypeElement,U,posgp',NameFileMesh,DOFl');
    [Q,R]=qr(SNAP_cluster.u.U);
    A(:,:,i)=Q(:,:);
    i=i+1;
end
C=A(:,:,1).'*A(:,:,2);
[Q,W,E]=svds(C,8);
disp(diag(W));
