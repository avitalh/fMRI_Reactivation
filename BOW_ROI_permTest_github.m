% BOW ROI sign-flip permutation test
clear all

ISC=1; % 1 for between-participant analysis, 0 for within-participant analysis

% set number of diagonals to remove by analysis type
if ISC~=1
    numDiag=10;
else
    numDiag=1;
end

% load BOW similarity matrix
load('BOW_similarity_matrix');

% load minimum word matrix
load('min_words');

% regress the minWord matrix from the BOW matrix
X = [ones(length(LangCorrs)*length(LangCorrs),1) minWord(:)];
[b,bint,res] = regress(LangCorrs(:),X);
LangCorrs=reshape(res,length(LangCorrs),length(LangCorrs));
clear b bint idx X

%remove diagonals
m=length(LangCorrs);
if numDiag>1
    X = full(spdiags(bsxfun(@times,ones(m,1),nan(1,numDiag)),[(numDiag/2-1)*-1:numDiag/2],m,m));
else
    X = eye(m);
    X(X==1)=nan;
end

LangCorrs=LangCorrs+X;

noRep=[2,18,21,42,27]; %short scenes with no representation (<=10 TRs), 27 is the last scene of the first scan - no representation
cartoon=[1,28]; %1/28 is cartoon
short=[3:15]; %scenes that are too short to have a reliable BOW inter-scene similarity
removeScenes=sort(unique([noRep,cartoon,short])); % short scenes and those with no BOW rep or no EB

% make matrix symmetrical
LangCorrs(:,removeScenes)=nan;
LangCorrs(removeScenes,:)=nan;

% load ROI representations (e.g. for the sherlock dataset in the between-participant 
% analysis, load the 21st year PCUN ROI defined from the between-participant analysis)
load('sherlockData_21_PCUN_ISC_representations');

for iSub=1:17
    if ISC~=1
        boundScene=corr(squeeze(repEB(:,:,iSub)),squeeze(sceneAvg(:,:,iSub)),'rows','pairwise');
    else
        ss=squeeze(sceneAvg(:,:,iSub));
        group=repEB;
        group(:,:,iSub)=[];
        boundScene=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
    end
    
    boundScene=boundScene+X;
    boundScene(:,removeScenes)=nan; %for past>future symmetry
    boundScene(removeScenes,:)=nan; %for past>future symmetry
    
    Past(iSub)=dist_and_fisher(getTriangilar(boundScene,1),getTriangilar(LangCorrs,1),'correlation');
    Future(iSub)=dist_and_fisher(getTriangilar(boundScene,0),getTriangilar(LangCorrs,0),'correlation');
end

% reactivation index
reactivationIdx=Past-Future; 

% permutation test
for itr=1:1000
    idx=randi([0, 1], 1,length(reactivationIdx)); %sign-flip index
    idx(idx==0)=-1;
    Rand(:,itr)=reactivationIdx.*idx;
end

% the test statistic
[h p c t]=ttest(reactivationIdx);
testStatistic=t.tstat;

% random statistics
[h p c t]=ttest(Rand);
t=t.tstat;

% calculate one-tailed p-value
pval=1-sum(testStatistic>=t)/length(t);
if pval==0
    pval=1/(length(t)+1);
end
