clear;

addpath("Postprocess\")

%%Postproces Modes GID
files=dir("GIDPOST/*NSS*");
DATA=0;

for file=files'
    if(isfile(strcat("GIDPOST\",file.name,'/',file.name,'_MODES.res')) )
        disp("Already computed modes");
    else
        load(strcat("GIDPOST\",file.name,'/',file.name,".mat"));
        U=SNAP_cluster.u.U(DOFl,:);
        GidPostProcessModes(COOR,CN,TypeElement,U,posgp',NameFileMesh,DOFl');
    end  
end
%%
%%Compute all subspace calculations for all possible permutations
%As long as subspace(A)>=subspace(B)
disp("Calculating angles between subspaces LidDriven")
files=dir("GIDPOST/LidDriven*NSS*");
combinations=nchoosek(files,2);
for i=1:size(combinations,1)
    load(strcat("GIDPOST\",combinations(i,1).name,'/',combinations(i,1).name,".mat"));
    U1=SNAP_cluster.u.U;
    load(strcat("GIDPOST\",combinations(i,2).name,'/',combinations(i,2).name,".mat"));
    U2=SNAP_cluster.u.U;
    if(size(U1,2)>=size(U2,2)) %p>=q
        [Q1,R1]=qr(U1,"econ");
        [Q2,R2]=qr(U2,"econ");
        C=Q1.'*Q2;
        [Y,SIGMA,Z]=svd(C);
        cosangles=diag(SIGMA);
        disp(strcat("Angles between ",combinations(i,1).name," and ",combinations(i,2).name))
        disp(cosangles);
    else     
        disp("Not valid dimensions");
    end
end
%%
disp("Calculating angles between subspaces Cylinder")
files=dir("GIDPOST/Cylinder40NSS*");
combinations=nchoosek(files,2);
for i=1:size(combinations,1)
    load(strcat("GIDPOST\",combinations(i,1).name,'/',combinations(i,1).name,".mat"));
    U1=SNAP_cluster.u.U;
    load(strcat("GIDPOST\",combinations(i,2).name,'/',combinations(i,2).name,".mat"));
    U2=SNAP_cluster.u.U;
    if(size(U1,2)>=size(U2,2)) %p>=q
        [Q1,R1]=qr(U1,"econ");
        [Q2,R2]=qr(U2,"econ");
        C=Q1.'*Q2;
        [Y,SIGMA,Z]=svd(C);
        cosangles=diag(SIGMA);
        disp(strcat("Angles between ",combinations(i,1).name," and ",combinations(i,2).name))
        disp(cosangles);
    else     
        disp("Not valid dimensions");
    end
end