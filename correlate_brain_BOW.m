function correlate_brain_BOW

ISC=1; % 1 for between-participant analysis, 0 for within-participant analysis

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

% load BOW similarity matrix
load('BOW_similarity_matrix');

% load minimum word matrix
load('min_words');

short=[1,2,18,21,42]; % short scenes (<10 TRs), no valid representation
irrelevant=[1,27,28]; % scenes 1,28 are cartoons, scene 27 is the end of the first scan and has no valid event boundary representation within the data
removeScenes=unique([short,irrelevant]);

% remove these scenes from the BOW and minimum word matrices
LangCorrs(:,removeScenes)=nan;
LangCorrs(removeScenes,:)=nan;
minWord(:,removeScenes)=nan;
minWord(removeScenes,:)=nan;

% remove diagonals from the BOW/minWord matrices
m=size(LangCorrs,1);
if numDiag>1
    X = full(spdiags(bsxfun(@times,ones(m,1),nan(1,numDiag)),[(numDiag/2-1)*-1:numDiag/2],m,m)); %create matrix of zeros with nans at the diagonals
else
    X = eye(m); 
    X(X==1)=nan;
end

LangCorrs=LangCorrs+X;
minWord=minWord+X;

% get representation file names
files=dir('*_reps.mat');

% load and concatenate representations across subjects
for iSub=1:length(files)
    load(files(iSub).name,'repEB','sceneAvg','sphereCenters');
    DataEB(:,:,:,iSub)=repEB;
    DataScenes(:,:,:,iSub)=sceneAvg;
end
clear repEB sceneAvg

% now create single-subject maps
for iSub=1:size(DataEB,3) %length(files)
    [mapT1] = deal(zeros(size(template.img))); % preallocate map files
    
    for iSphere=1:length(sphereCenters.x)
        % get current sphere data
        repEB=squeeze(DataEB(:,:,iSphere,:));
        sceneAvg=squeeze(DataScenes(:,:,iSphere,:));
        
        %remove padding (marked by 999 values) in small spheres
        temp=squeeze(repEB(:,:,iSub));
        idx=find(nanmean(temp,2)==999);
        repEB(idx,:,:)=[];
        sceneAvg(idx,:,:)=[];
        
        % remove unreliable matrices (>66.6% missing values in either part of the scan)
        temp=squeeze(repEB(:,:,iSub));
        if sum(~isnan(nanmean(temp(:,1:27)')))<41 || sum(~isnan(nanmean(temp(:,28:end)')))<41
            mapT1(sphereCenters.x(iSphere),sphereCenters.y(iSphere),sphereCenters.z(iSphere))=nan;
            continue
        end
        
        if ISC~=1 % within-participant analysis
            boundSceneE=corr(squeeze(repEB(:,:,iSub)),squeeze(sceneAvg(:,:,iSub)),'rows','pairwise');
            
        else % between-participant analysis
            ss=squeeze(sceneAvg(:,:,iSub));
            group=repEB;
            group(:,:,iSub)=[];
            boundSceneE=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
        end
        
        % remove diagonals for the neural similarity matrix
        boundSceneE=boundSceneE+X;
        
        % To make the matrix symmetrical, we must remove both rows and columns of removed representations
        boundSceneE(removeScenes,:)=nan;
        boundSceneE(:,removeScenes)=nan;
        
        % corrlate brain-BOW matrices and fisher transform
        PastEB=dist_and_fisher(getTriangular(boundSceneE,1),getTriangular(LangCorrs,1),'correlation',getTriangular(minWord,1));
        FutureEB=dist_and_fisher(getTriangular(boundSceneE,0),getTriangular(LangCorrs,0),'correlation',getTriangular(minWord,0));
        
        % claculate reactivation index and assign to map
        mapT1(sphereCenters.x(iSphere),sphereCenters.y(iSphere),sphereCenters.z(iSphere))=PastEB-FutureEB;
        
    end
end
