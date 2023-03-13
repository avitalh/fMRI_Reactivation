path='dataPath/'; %placeholder
pathO='outputPath/'; %placeholder
ROIpath='ROIpath/';

ISC=0; % 1 for between-participant analysis, 0 for within-participant analysis
iLag=1; % number of TRs needed to shift to account for the HRF delay (the Sherlock dataset already contains a shift of 3 TRs, so we're adding one more TR)

if ISC==1
    numDiag=1;
else
    numDiag=10;
end

% load a list of event boundaries 
load([path 'allBoundaries'],'allBoundaries');
DiffBoundaries=[allBoundaries(1);DiffBoundaries]; % length of scenes
sceneBegin=[1;allBoundaries(1:end-1)+1]; % beginnings of scenes

% find functional data files
files=dir([path '*.nii']);

minSceneLen=20; % minimal number of TRs in each scene
preEB=8; % number of TRs to plot before the event boundary (4 TRs with regards to the boundaries uncorrected for HRF delay)
postEB=8; % number of TRs to plot after the event boundary (12 TRs with regards to the boundaries uncorrected for HRF delay)

% load ROI mask
ROI=load_nii([ROIpath '21_Pcun_ROI_ISC.nii']); % for the Sherlock dataset, use 21st year ROI
ROI=ROI.img;

% get the coordinates of ROI voxels
idx=find(ROI>0);
[i j k]=ind2sub(size(ROI),idx);
IDX=[i j k];
clear i j k idx

% preallocate variables 
sceneRep=nan(length(files),length(allBoundaries)-1,size(IDX,1),preEB+postEB); % representations of single TRs around the boundaries
sceneAvg=nan(length(files),length(allBoundaries)-1,size(IDX,1)); % avg. scene representations 

% get representations of scenes and TRs
for iSub=1:length(files)
    clear tcs corrs data 
    
    % load functional data
    nii=load_nii([path files(iSub).name]);
    data=nii.img; 
    clear nii
    
    % extract ROI timecourses
    tcs=zeros(length(IDX),size(data,4));
    for j=1:length(IDX)
        tcs(j,:)=data(IDX(j,1),IDX(j,2),IDX(j,3),:);
    end
    
    % remove voxels with no signal
    idx=find(mean(tcs')==0);
    tcs(idx,:)=[];
    
    for iScene=1:length(allBoundaries)-1 % not enouugh TRs after last event boundary
        if DiffBoundaries(iScene)>=minSceneLen
            sceneRep(iSub,iScene,:,:)=tcs(:,allBoundaries(iScene)-preEB+1+iLag:allBoundaries(iScene)+postEB+iLag); % TR representations
            sceneAvg(iSub,iScene,:)=mean(tcs(:,sceneBegin(iScene)+5+iLag:allBoundaries(iScene)-5+iLag),2); % averaged scene represenattion
        else
            sceneRep(iSub,iScene,:,:)=nan;
            sceneAvg(iSub,iScene,:)=nan;
        end
    end
end

% correlate the single-TR representations with the averaged scene representations
for iSub=1:length(files)
    for iEvent=1:length(allBoundaries)-1
        for jEvent=1:length(allBoundaries)-1
            if ISC~=1
                Corr(iEvent,jEvent,:)=corr(squeeze(sceneRep(iSub,iEvent,:,:)),squeeze(sceneAvg(iSub,jEvent,:)),'rows','pairwise');
            else
                ss=squeeze(sceneAvg(iSub,jEvent,:));
                group=squeeze(sceneRep(:,iEvent,:,:));
                group(iSub,:,:,:)=[];
                Corr(iEvent,jEvent,:)=corr(squeeze(nanmean(group,1)),ss,'rows','pairwise');
            end
        end
    end
    
    %remove diagonals
    m=size(Corr,1);
    if numDiag>1
        X = full(spdiags(bsxfun(@times,ones(m,1),nan(1,numDiag)),[(numDiag/2-1)*-1:numDiag/2],m,m)); %create matrix of zeros with nans at the diagonals
    else
        X = eye(m);
        X(X==1)=nan;
    end
    X=repmat(X,1,1,size(Corr,3));
    Corr=Corr+X;
    
    % get past correlations
    idx=tril(ones(length(allBoundaries)-1,length(allBoundaries)-1));
    idx(logical(eye(length(allBoundaries)-1)))=0;
    Past=Corr;
    Past=Past.*repmat(idx,1,1,preEB+postEB);
    Past(Past==0)=nan;
    
    % get the mean correlation with all past scenes for each TR position
    % relative to the boundary (e.g. the TR immediatly preceding the boundary)
    for i=1:preEB+postEB
        temp=squeeze(Past(:,:,i));
        Past_all(iSub,i)=mean(temp(~isnan(temp(:))));
    end
    
    
    % get future correlations
    idx=triu(ones(length(allBoundaries)-1,length(allBoundaries)-1));
    idx(logical(eye(length(allBoundaries)-1)))=0;
    Future=Corr;
    Future=Future.*repmat(idx,1,1,preEB+postEB);
    Future(Future==0)=nan;
    
    % get the mean correlation with all future scenes for each TR position
    % relative to the boundary (e.g. the TR immediatly preceding the boundary)
    for i=1:preEB+postEB
        temp=squeeze(Future(:,:,i));
        Future_all(iSub,i)=mean(temp(~isnan(temp(:))));
    end
    
end

% get the reactivation index
Diff=Past_all-Future_all;

%% plot

figure(1); clf
axis('equal')
axis square
H=shadedErrorBar([],mean(Diff,1),std(Diff,1)./sqrt(size(Past_all,1)-1),'k',1)
xlabel('within-scene TRs')
ylabel('Reactivation index')
hold on
vline(4,'r') % 4 TRS before EB (scene TRs) not accounting for HRF delay
hline(0,'k');
xlim([1,16])
set(gca,'Xtick',0:4:16,'XTickLabel',-4:4:12)

%% statistical analysis

% get single-subject effect sizes
for iSub=1:size(Diff,1)
    signal=Diff(iSub,preEB-4+1:preEB-4+6); % TRs afetr uncorrected event boundaries (9sec)
    baseline=Diff(iSub,preEB-4-1:preEB-4);% TRs prior to uncorrected event boundaries 
    
    ssVals(iSub)=sum(signal)-sum(baseline); 
end

% sign flip permutation test (two-sided)
[h p c t]=ttest(ssVals);
testStatistic=t.tstat;

for itr=1:1000
    idx=randi([0, 1], 1,length(ssVals)); %sign-flip index
    idx(idx==0)=-1;
    Rand(:,itr)=ssVals.*idx;
end
    
% random statistics
[h p c t]=ttest(Rand);
t=t.tstat;

% calculate p-value
pval=1-sum(testStatistic>=t)/length(t);
if pval==0
    pval=1/(length(t)+1);
end

% two-tailed correction
if pval<=0.05
    pval=pval*2;
end 
