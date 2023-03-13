% create data for intra-subject analysis
clear
addpath spm12\external\fieldtrip\preproc % SPM directory with a high-pass filter function (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

%simulation parameters
TR=1500;
fs=1/(TR/1000);
filter_width=140; % as in the preprocessed Sherlock data
N=1976;  % number of datapoints
Lag=1; % number of TRs needed to shift to account for the HRF delay (the Sherlock dataset already contains a shift of 3 TRs, so we're adding one more TR)
numDiag=10; %number of diagonals to remove from the analysis
distCtr=10; % distance of control timepoints from event boundaries
sceneDist=5; % distance of scene representations from adjecent event boundaries
numSub=10; % number of simulated subjects
numItr=1000; %number of iterations

load('allBoundaries'); % load a variable called allBoundaries, which is a list of event boundaries (change path)
sceneBegin=[1;allBoundaries(1:end-1)+1]; % scene beginning points
DiffBoundaries=diff(allBoundaries); % scene lengths

% preallocate variables
boundSceneE=nan(length(allBoundaries)-1,length(allBoundaries)-1,numSub,numItr); % event boundary correlations
boundSceneP=nan(length(allBoundaries)-1,length(allBoundaries)-1,numSub,numItr); % past control correlations
boundSceneF=nan(length(allBoundaries)-1,length(allBoundaries)-1,numSub,numItr); % future control correlations

for itr=1:numItr
    for iSub=1:numSub % index of subjects
        clear data
        
        % simulate 1/f signals
        for i=1:120 %index of voxels
            y = pinknoise(N);
            data(i,:)=y;
        end
        
        % filter the data
        [data, B, A] = ft_preproc_highpassfilter(data,fs,1/filter_width); % SPM function 
                
        %get representations
        for iEvent=1:length(allBoundaries)-1 % all events except the last one (as it has no future control timepoint within the dataset)
            eventBound=allBoundaries(iEvent)+Lag; % index of current event shifted to account for the HRF delay
            repE(:,iEvent)=data(:,eventBound); % representation at event boundary
            repP(:,iEvent)=data(:,eventBound-distCtr); % representation at past control timepoint
            repF(:,iEvent)=data(:,eventBound+distCtr); % representation at future control timepoint
        end
        
        % get mean scene representataions (shifted by sceneDist and iLag)       
        for iEvent=1:length(allBoundaries)-1
            repScene(:,iEvent)=mean(data(:,sceneBegin(iEvent)+sceneDist+Lag:allBoundaries(iEvent)-sceneDist+Lag),2); 
        end
        
        % calculate correlation matrices
        boundSceneE(:,:,iSub,itr)=corr(squeeze(repE),squeeze(repScene),'rows','pairwise');
        boundSceneP(:,:,iSub,itr)=corr(squeeze(repP),squeeze(repScene),'rows','pairwise');
        boundSceneF(:,:,iSub,itr)=corr(squeeze(repF),squeeze(repScene),'rows','pairwise');
    end
end


%% plot
close all
figure(1); clf

% remove bad scenes for permutaton (according to the true data)
short=[1,2,18,21,42, 27]; % short scenes, no valid representation, and scene 27 is the end of the first scan and has no valid event boundary representation within the data
cartoon=[1,28]; % scenes 1,28 are cartoons
removeScenes=unique([short,cartoon]);

boundSceneE(removeScenes,:,:,:)=nan;
boundSceneE(:,removeScenes,:,:)=nan;
boundSceneP(removeScenes,:,:,:)=nan;
boundSceneP(:,removeScenes,:,:)=nan;
boundSceneF(removeScenes,:,:,:)=nan;
boundSceneF(:,removeScenes,:,:)=nan;

% plot average reactivation index matrix (across iterations and subjects)
EB=squeeze(mean(mean(boundSceneE,3),4));

%remove diagonals (X masks diagonals)
m=size(boundSceneE,1);
X = full(spdiags(bsxfun(@times,ones(m,1),nan(1,numDiag)),[(numDiag/2-1)*-1:numDiag/2],m,m)); %create matrix of zeros with nans at the diagonals
EB_d=EB+X;

% plot the full matrix to visualize the autocorrelation bias
subplot(2,2,1) 
imagesc(EB)
axis square
title('full matrix')
caxis([-0.1 0.1])
colorbar
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)

% plot the matrix without the impacted diagonals to visualize the removal of the bias
subplot(2,2,2) 
imagesc(EB_d)
axis square
title('matrix without diagonals')
caxis([-0.1 0.1])
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)

% create a reactivation index distributions to make sure the bias is indeed removed
for itr=1:size(boundSceneE,4) % iteration index
    for iSub=1:size(boundSceneE,3) % subject index
        temp=squeeze(boundSceneE(:,:,iSub,itr));
        idxE(iSub)=mean(getTriangilar(temp,1))-mean(getTriangilar(temp,0)); % reactivation index for event boundaries
        
        temp=squeeze(boundSceneP(:,:,iSub,itr));
        idxP(iSub)=mean(getTriangilar(temp,1))-mean(getTriangilar(temp,0)); % reactivation index for past control point
        
        temp=squeeze(boundSceneF(:,:,iSub,itr));
        idxF(iSub)=mean(getTriangilar(temp,1))-mean(getTriangilar(temp,0)); % reactivation index for future control point
    end
    
    % average across subjects to geta group index
    groupE(itr)=mean(idxE);
    groupP(itr)=mean(idxP);
    groupF(itr)=mean(idxF);
end

subplot(2,2,3) 
histogram(groupE,100,'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
xline(mean(groupE),':r','linewidth',2.5)
title('reactivation index, full matrix')
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)
box off
axis square

%remove diagonals
boundSceneE=boundSceneE+X;
boundSceneP=boundSceneP+X;
boundSceneF=boundSceneF+X;

% create the new distribution
for itr=1:size(boundSceneE,4)
    for iSub=1:size(boundSceneE,3)
        temp=squeeze(boundSceneE(:,:,iSub,itr));
        idxE(iSub)=mean(getTriangilar(temp,1))-mean(getTriangilar(temp,0));
        
        temp=squeeze(boundSceneP(:,:,iSub,itr));
        idxP(iSub)=mean(getTriangilar(temp,1))-mean(getTriangilar(temp,0));
        
        temp=squeeze(boundSceneF(:,:,iSub,itr));
        idxF(iSub)=mean(getTriangilar(temp,1))-mean(getTriangilar(temp,0));
    end
    
    groupE(itr)=mean(idxE);
    groupP(itr)=mean(idxP);
    groupF(itr)=mean(idxF);
end

subplot(2,2,4) 
histogram(groupE,100,'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
xline(mean(groupE),':r','linewidth',2.5)
title('reactivation index, without diagonals')
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)
box off
axis square
set(gcf,'color','w');

% now plot the histograms of control points (after removal of diagonals), and of the balanced control
figure(2); clf
subplot(1,3,1)
histogram(groupP,100,'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
xline(mean(groupP),':r','linewidth',2.5)
title('past control')
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)
box off
axis square
set(gcf,'color','w');

subplot(1,3,2)
histogram(groupF,100,'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
xline(mean(groupF),':r','linewidth',2.5)
title('future control')
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)
box off
axis square
set(gcf,'color','w');

subplot(1,3,3)
balancedControl=groupE-(groupP+groupF)/2;
histogram(balancedControl,100,'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
xline(mean(balancedControl),':r','linewidth',2.5)
title('balanced control')
set(gca,'linewidth',2.5)
set(gca,'fontsize',15)
box off
axis square
set(gcf,'color','w');