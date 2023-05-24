function [inELEM,onELEM,TriLocal] = IsInsideGeneral2023(xLOC,COOR,CN,elemLOC,IND_POLYG,POLYINFO)
TriLocal = [] ; 
ndim = size(COOR,2) ;
CNloc = CN(elemLOC,:); % Nodes forming the element
if ndim == 2
    INDnodes =CNloc(IND_POLYG) ;  % Corner nodes (to define a polygon)
    COORelem = COOR(INDnodes,:) ; % Coordinates of the elemen
    [inELEM,onELEM] = inpolygon(xLOC(1),xLOC(2),COORelem(:,1),COORelem(:,2)) ;
else
    % Use triangulation
    % ------------------
    INDnodes =CNloc(IND_POLYG) ;  % Corner nodes
    INDnodes = INDnodes(1:end-1) ;
    COORelem = COOR(INDnodes,:) ; % Coordinates of the elemen
  %  error('REvise this ! ')
    if isempty(POLYINFO.TriangulationDelaunay{elemLOC})
    TriLocal = delaunayTriangulation(COORelem);
    else
        TriLocal = POLYINFO.TriangulationDelaunay{elemLOC} ; 
    end
    CHECKINSIDE = pointLocation(TriLocal,xLOC) ;
    if isnan(CHECKINSIDE)
        inELEM = 0 ; onELEM = 0 ;
    else
        inELEM = 1; onELEM = 0 ;
    end    
end