function [wST,wSTs] = DetermineWeightsST(detJeALL,weigREP,ngaus,g,nelem,wST,wSTs,nstrain,ROWSglo);

wLOCa = detJeALL.*weigREP(g:ngaus:nelem*ngaus) ;
wSTs(g:ngaus:nelem*ngaus) =  wLOCa ;
%wLOCa = [5 6]' ;
wLOCb = repmat(wLOCa',nstrain*2,1) ;
wLOCb  = wLOCb(:) ;
wST(ROWSglo) = wLOCb(:);

end