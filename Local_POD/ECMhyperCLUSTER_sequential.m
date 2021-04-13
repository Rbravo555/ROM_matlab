function [ECMdata_cluster,setCandidates] = ECMhyperCLUSTER_sequential(BasisU_cluster,BasisPone_cluster,BasisStwo_cluster,...
    OPERFE,DISP_CONDITIONS,DATA,DATAoffline)

if nargin  == 0
    load('tmp.mat')
end


nclusters = length(BasisU_cluster) ;

ECMdata_cluster = cell(1,nclusters) ;

disp('SEQUENTIAL METHOD FOR HYPERREDUCTION')
disp('********************************************************************')
disp(['Determining basis matrices for internal forces for each cluster'])
disp('********************************************************************')
BasisFint_cluster = cell(size(BasisStwo_cluster)); 
DOFl = DISP_CONDITIONS.DOFl ;
DOFr = DISP_CONDITIONS.DOFr ;
for icluster = 1:nclusters
    
    BstRED_l = OPERFE.Bst(:,DOFl)*BasisU_cluster{icluster} ; 
    BasisFint_cluster{icluster} = QbasisMatrixIntegrand(BstRED_l,BasisPone_cluster{icluster},DATA,OPERFE.wSTs,DATAoffline) ; 
   disp(['CLUSTER = ',num2str(icluster),';  N modes FINT =  ',num2str(size(BasisFint_cluster{icluster},2))]) ;   
    
    %     BasisPone = BasisPone_cluster{icluster} ;
    %     BasisStwo = BasisStwo_cluster{icluster} ;
    %     ECMdata_cluster{icluster} = ECM_clustersNC(OPERFE,DISP_CONDITIONS,BasisU,DATAoffline,BasisPone,DATA,BasisStwo) ;
    %     disp(['*****************************'])
end

% Sort the number of columns of each BasisFint_cluster 
DATAoffline = DefaultField(DATAoffline,'OrderSequentialECM','SORTED') ; % If =2, 
[dummy,nmodesFINT] = cellfun(@size,BasisFint_cluster) ; 

switch DATAoffline.OrderSequentialECM  
    case 'SORTED'
[III,NEW_ORDER_clusters] = sort(nmodesFINT,'descend') ; 
    case 'RANDOM'
    NEW_ORDER_clusters = randperm(length(BasisFint_cluster)) ; 
    case 'RANDN_fixed'
        rng(1) 
         NEW_ORDER_clusters = randperm(length(BasisFint_cluster)) ; 
  %  NEW_ORDER_clusters = 1:length(BasisFint_cluster) ; 
end

% Next we run a loop over all clusters, starting with the cluster with
% higher number of modes (i.e., points)
disp('*************************************************++')
disp(['Sequential ECM ...'])
setCandidates = [] ; 
for iclusterLOC = 1:nclusters
    
   
    icluster = NEW_ORDER_clusters(iclusterLOC) ; 
     if  icluster == 8
         disp('borrar esto')
    end
    
    disp('*************************************************++')
    disp(['CLUSTER = ',num2str(icluster)])
    disp('*************************************************++')
    BasisStwo = BasisStwo_cluster{icluster} ;
    [ECMdata_cluster{icluster},setCandidates ]= ECM_clustersNC(OPERFE.wSTs,DATAoffline,BasisFint_cluster{icluster},DATA,BasisStwo,setCandidates) ;
    disp(['*****************************'])
    
    
    disp(['*********************************************************************'])
disp(['Number of displacement modes = ',num2str(size(BasisU_cluster{icluster},2))])
disp(['Number of PK2-stress modes = ',num2str(size(BasisStwo_cluster{icluster},2))])
disp(['Number of PK1-stress modes = ',num2str(size(BasisPone_cluster{icluster},2))])
end

disp(['Total number of points (had they been  selected independently) = ',num2str(sum(nmodesFINT))])

disp(['Total number of ECM candidate points = ',num2str(length(setCandidates))])


setElements = large2smallREP(setCandidates,DATA.MESH.ngaus) ;
disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements (for all clusters)'])

 disp(num2str(setElements'))
    clipboard('copy',num2str(setElements'));