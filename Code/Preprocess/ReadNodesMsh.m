function [rnod] = ReadNodesMsh(NameFileMesh,ILINE)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
NAMEFILE=NameFileMesh(1:end-4);
MeshName=split(NAMEFILE,'/',1);
nameDAT3 = [NAMEFILE,'.gid/',char(MeshName(end)),'-2.dat'] ;
NODES_Node = {} ;
if exist(nameDAT3,'file')
    DDD = load(nameDAT3) ;
    if ~isempty(DDD)
        LABELS_Nods = unique(DDD(:,2)) ;
        nmaxLABEL = max(LABELS_Nods) ;
        NODES_Node = cell(1,nmaxLABEL) ;
        for inodLOC = 1:length(LABELS_Nods)
            inod = LABELS_Nods(inodLOC) ;
            INDX = find(DDD(:,2) == inod) ;
            NODES_Nods{inodLOC} = DDD(INDX,1) ;
        end
        rnod=cell2mat(NODES_Nods);
    else
        rnod=[];
    end
    
end



end