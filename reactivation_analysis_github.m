% Creates a reactivation index whole-brain map
function reactivation_analysis
pathO='outputPath'; %placeholder

ISC=0; % 1 for between-participant analysis, 0 for within-participant analysis

% remove diagonals by analysis type
if ISC==1
    numDiag=1;
    suffix='ISC';
else
    numDiag=10;
    suffix='intra';
end

%load template MNI map
template=load_nii('template_map.nii'); %load template map

short=[1,2,18,21,42]; % short scenes (<10 TRs), no valid representation
irrelevant=[1,27,28]; % scenes 1,28 are cartoons, scene 27 is the end of the first scan and has no valid event boundary representation within the data
removeScenes=unique([short,irrelevant]);
clear short irrelevant

files=dir('*reps.mat'); % get representation file names

% first load and concatenate representations across subjects
for iSub=1:length(files)
    % load representations
    load(files(iSub).name,'repEB','sceneAvg','repP','repF','sphereCenters');
    
    DataEB(:,:,:,iSub)=repEB;
    DataScenes(:,:,:,iSub)=sceneAvg;
    DataP(:,:,:,iSub)=repP;
    DataF(:,:,:,iSub)=repF;
end
clear repEB sceneAvg repP repF

% now create single-subject maps
for iSub=1:length(files)
    [mapT1,mapT2] = deal(zeros(size(template.img))); % preallocate map files
    
    for iSphere=1:length(sphereCenters.x)
        % get current sphere data
        repEB=squeeze(DataEB(:,:,iSphere,:));
        sceneAvg=squeeze(DataScenes(:,:,iSphere,:));
        repP=squeeze(DataP(:,:,iSphere,:));
        repF=squeeze(DataF(:,:,iSphere,:));
        
        %remove padding (marked by 999 values) in small spheres
        temp=squeeze(repEB(:,:,iSub));
        idx=find(nanmean(temp,2)==999);
        repEB(idx,:,:)=[];
        sceneAvg(idx,:,:)=[];
        repP(idx,:,:)=[];
        repF(idx,:,:)=[];
        
        % ignore unreliable sphere representations (spheres < 45 voxels)
        if size(repEB,1)<45
            continue
        end
        
        if ISC~=1 % within-participant analysis
            boundSceneE=corr(squeeze(repEB(:,:,iSub)),squeeze(sceneAvg(:,:,iSub)),'rows','pairwise');
            boundSceneP=corr(squeeze(repP(:,:,iSub)),squeeze(sceneAvg(:,:,iSub)),'rows','pairwise');
            boundSceneF=corr(squeeze(repF(:,:,iSub)),squeeze(sceneAvg(:,:,iSub)),'rows','pairwise');
            
        else % between-participant analysis
            ss=squeeze(sceneAvg(:,:,iSub)); % single-subject scene representations
            group=repEB; % group boundary representations
            group(:,:,iSub)=[]; % exclude the subject whose scene representations will be used
            boundSceneE=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
            
            group=repP;
            group(:,:,iSub)=[];
            boundSceneP=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
            
            group=repF;
            group(:,:,iSub)=[];
            boundSceneF=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
        end
        
        % remove diagonals
        m=size(boundSceneE,1);
        if numDiag>1
            X = full(spdiags(bsxfun(@times,ones(m,1),nan(1,numDiag)),[(numDiag/2-1)*-1:numDiag/2],m,m)); %create matrix of zeros with nans at the diagonals
        else
            X = eye(m);
            X(X==1)=nan;
        end
        
        boundSceneE=boundSceneE+X;
        boundSceneP=boundSceneP+X;
        boundSceneF=boundSceneF+X;
        
        % To make the matrix summetrical, we must remove the rows (and not only the columns) of invalid representations
        boundSceneE(removeScenes,:)=nan;
        boundSceneP(removeScenes,:)=nan;
        boundSceneF(removeScenes,:)=nan;
        
        % Fisher correction
        boundSceneE=0.5*(log(1+boundSceneE)-log(1-boundSceneE));
        boundSceneP=0.5*(log(1+boundSceneP)-log(1-boundSceneP));
        boundSceneF=0.5*(log(1+boundSceneF)-log(1-boundSceneF));
        
        %calculate the reactivation index
        indexPastFutureE=nanmean(getTriangular(boundSceneE,1))-nanmean(getTriangular(boundSceneE,0));
        indexPastFutureP=nanmean(getTriangular(boundSceneP,1))-nanmean(getTriangular(boundSceneP,0));
        indexPastFuturF=nanmean(getTriangular(boundSceneF,1))-nanmean(getTriangular(boundSceneF,0));
        
        % assign index to map
        mapT1(sphereCenters.x(iSphere),sphereCenters.y(iSphere),sphereCenters.z(iSphere))=indexPastFutureE;
        mapT2(sphereCenters.x(iSphere),sphereCenters.y(iSphere),sphereCenters.z(iSphere))=indexPastFutureE-(indexPastFutureP+indexPastFuturF)/2;
        
    end
    
    save([pathO files(iSub).name(1:end-4) '_' suffix '_reactivation'],'mapT1','mapT2');
end