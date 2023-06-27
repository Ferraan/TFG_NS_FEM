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

