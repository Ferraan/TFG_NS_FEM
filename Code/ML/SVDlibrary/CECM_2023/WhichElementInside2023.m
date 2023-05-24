function [elemCONTAINER,POLYINFO] = WhichElementInside2023(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG,POLYINFO,inew) ;
if nargin == 0
    load('tmp.mat')
end

ielem = 1;
elemCONTAINER = [] ;
[ELEMnear, aaaa ]= find(VAR_SMOOTH_FE.CN == INDnear) ;  % Elements sharing INDnear.

if ~isempty(POLYINFO.setElements)
    currentELEMENT = POLYINFO.setElements(inew) ;
    
    if currentELEMENT ~=0
        % This is the last element stored in memory. We shall search first in
        % this element, as well as on the neighboring elements
        nNEIGHS = VAR_SMOOTH_FE.CONNECT_info.ElemShared(currentELEMENT)  ;
        NEIGH_elemes = VAR_SMOOTH_FE.CONNECT_info.TableElements(currentELEMENT,1:nNEIGHS) ;
        ELEMnear = setdiff(ELEMnear,currentELEMENT) ;
        ELEMnear = unique([ELEMnear;NEIGH_elemes(:)],'stable') ;
        % Finally
        ELEMnear = [currentELEMENT;ELEMnear] ;
        
    end
    
end

while ielem <= length(ELEMnear)
    elemLOC = ELEMnear(ielem) ;
    
    [inELEM,onELEM,TriLocal] = IsInsideALL2023(xLOC,VAR_SMOOTH_FE,elemLOC,IND_POLYG,POLYINFO)  ;
    POLYINFO.TriangulationDelaunay{elemLOC} = TriLocal ;
    if inELEM == 1 || onELEM == 1
        elemCONTAINER  = elemLOC ;
        %             if elemCONTAINER == 99
        %                 disp('Borrar esto')
        %             end
        break
    end
    ielem = ielem + 1;
end

% if isempty(elemCONTAINER)
%     disp('')
% end