function [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED,...
    COOR,CN,posgp,TypeElement,DOFl,NameFileMesh)

if nargin == 0
    load('tmp.mat')
end

[ISCURRENT_CLUSTER,istepLOC] = ismember(istep,DATA.STORE.NSTEPS_CLUSTER{icluster}) ;

if ISCURRENT_CLUSTER
    
    %      if ~isfield(SNAP,'DISP')
    %          SNAP.DISP(:,istepLOC) = d;
    %      end
    %       if ~isfield(SNAP,'PK2STRESS')
    %          SNAP.PK2STRESS(:,istepLOC) = StwoST;
    %       end
    %       if ~isfield(SNAP,'GLSTRAINS')
    %          SNAP.GLSTRAINS(:,istepLOC) = EgreenlST;
    %       end
    
    if ~isstruct(CONVERGED) && CONVERGED == 1
        fff = fieldnames(SNAP) ;
        
        for iii = 1:length(fff)
            nrows = size(VAR.(fff{iii})) ;
            SNAP.(fff{iii})(1:nrows,istepLOC) = VAR.(fff{iii}) ;
        end
    end
    
    if istepLOC == length(DATA.STORE.NSTEPS_CLUSTER{icluster})  || CONVERGED == 0
        % STORE information previous cluster
        if DATA.STORE.COMPRESS_WITH_SVD  == 1
            disp(['Compressing snapshots cluster = ',num2str(icluster),'...'])
            abcd = tic ; 
            fff= fieldnames(SNAP) ;
            for iii = 1:length(fff)
                % [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA)
                DATASVD.RELATIVE_SVD = 1;
                
                [U,S,V] = RSVDT(SNAP.(fff{iii}),DATA.STORE.TOLERANCE_SVD_COMPRESSION.(fff{iii}),0,1,DATASVD) ;
                disp([fff{iii},' RANK = ',num2str(length(S))])
                
                %   SV = bsxfun(@times,V',S)';
                
                SNAP_cluster.(fff{iii}) = [] ;
                if ~isempty(S) > 0 
                SNAP_cluster.(fff{iii}).U = U ;
                SNAP_cluster.(fff{iii}).S = S ;
                SNAP_cluster.(fff{iii}).V = V ;
                else
                    % MOdification 10-APRIL-2023
                     SNAP_cluster.(fff{iii}).U = zeros(size(SNAP.(fff{iii}),1),1) ;
                SNAP_cluster.(fff{iii}).S = 0 ;
                SNAP_cluster.(fff{iii}).V = zeros(size(SNAP.(fff{iii}),2),1) ; 
                    
                end
            end
            TimeStore = toc(abcd) ;
            disp(['Done in ',num2str(TimeStore)])
            STEP_LOC = DATA.STORE.NSTEPS_CLUSTER{icluster}(1:istepLOC);
            save(DATA.STORE.NAME_MATFILE_STORE{icluster},'STEP_LOC','SNAP_cluster','COOR','CN','posgp','TypeElement','DOFl',"NameFileMesh") ;
        else
              STEP_LOC = DATA.STORE.NSTEPS_CLUSTER{icluster}(1:istepLOC);
            save(DATA.STORE.NAME_MATFILE_STORE{icluster},'STEP_LOC','SNAP') ;
         %   disp('Option not implemented')
        end
        % SAve information
        
      
    end
    
else
    
    fff= fieldnames(SNAP) ;
    % Change storing cluster
    icluster = icluster + 1;
    % Initialization new snapshot
    ntimestepsLOC = length(DATA.STORE.NSTEPS_CLUSTER{icluster}) ;
    for iii = 1:length(fff)
        ncomp = size(SNAP.(fff{iii}),1) ;
        SNAP.(fff{iii}) = zeros(ncomp,ntimestepsLOC) ;
    end
    
    [ISCURRENT_CLUSTER,istepLOC] = ismember(istep,DATA.STORE.NSTEPS_CLUSTER{icluster}) ;
    
    if CONVERGED == 1
        fff = fieldnames(SNAP) ;
        
        for iii = 1:length(fff)
            SNAP.(fff{iii})(:,istepLOC) = VAR.(fff{iii}) ;
        end
    end
end