% extracts rerpesentations
function getReps(iSub)
iLag=1; % how many TRs to shift to account for the HRF delay (the Sherlock dataset already contains a shift of 3 TRs, and we're adding one more TR)
numEvents=49; % the last event boundary has no representation (4 TR delay after end of scan)
distCtr=10; % distance of control timepoints from event boundaries

path='dataPath'; %placeholder
pathO='outputPath/'; %placeholder

% load MNI brain mask
nii=load_nii([path 'MNI_mask.nii']);
Mask=nii.img;

% load speher centers and sphere voxels
load([path 'spheres.mat']); % centers and sphere voxels

% load a list of event boundaries (allBoundaries)
load([path 'allBoundaries'],'allBoundaries');
DiffBoundaries=diff(allBoundaries); %compute scene lengths

%generate a list of control TRs
ctrP=allBoundaries-distCtr;
ctrF=allBoundaries+distCtr;

files=dir([path '*.nii']); % get functional file names

for iSub=1:length(files)
    nii=load_nii([path files(iSub).name]); % placeholder, get functional data of a single subject
    Data=single(nii.img);
    Data=(reshape(shiftdim(Data,3),size(Data,4),271633));
    
    % preallocate representation variables
    repEB=nan(123,numEvents,length(sphereCenters.x)); % representations at boundaries
    sceneAvg=nan(123,numEvents,length(sphereCenters.x)); % scene representations
    repP=nan(123,numEvents,length(sphereCenters.x)); % representations at past control points
    repF=nan(123,numEvents,length(sphereCenters.x)); % representations at future control points
    
    for iSphere=1:length(sphereInds)
        % sphere indices are relative to non-zero mask voxels (rather than to the whole MNI space),
        % so we need to convert them to MNI indices in order to extract relevant signals
        mask=zeros(size(Mask));
        nonZero=find(Mask);
        sphere=sphereInds{iSphere};
        mask(nonZero(sphere))=1;
        idx=find(~mask);
        data=Data';
        data(idx,:)=[];
        data=data';
        
        % extract the representations for boundaries and control points across the entire movie
        for iEvent=1:numEvents
            eventBound=allBoundaries(iEvent);
            temp=data(eventBound+iLag,:);
            % spheres can hold up to 123 voxels. Smaller spheres are padded with 999 to assign them to the representation matrix
            repEB(:,iEvent,iSphere)=padarray(temp,[0,123-length(temp)],999,'post');
            
            temp=data(ctrP(iEvent)+iLag,:); %past control
            repP(:,iEvent,iSphere)=padarray(temp,[0,123-length(temp)],999,'post');
            
            temp=data(ctrF(iEvent)+iLag,:); %future control
            repF(:,iEvent,iSphere)=padarray(temp,[0,123-length(temp)],999,'post');
        end
        
        repEB(:,[1,27,28],iSphere)=nan; %scenes 1/28 are cartoons, scene 27 ends the first scan - no event boundary representation in the data
        repP(:,[1,27,28],iSphere)=nan;
        repP(:,[1,27,28],iSphere)=nan;
        
        % get scene representations
        for iEvent=1:numEvents
            eventBound=allBoundaries(iEvent);
            if iEvent==1
                temp=mean(data(1+5+iLag:eventBound-5+iLag,:),1);
                sceneAvg(:,iEvent,iSphere)=padarray(temp,[0,123-length(temp)],999,'post');
            elseif DiffBoundaries(iEvent-1)>10 % 10 TRs is the minimum scene-length since we need to subtract 5TRs from the beginning and end of the scene
                temp=mean(data(allBoundaries(iEvent-1)+5+iLag:eventBound-5+iLag,:),1);
                sceneAvg(:,iEvent,iSphere)=padarray(temp,[0,123-length(temp)],999,'post');
            else
                sceneAvg(:,iEvent,iSphere)=nan;
            end
        end
        
        sceneAvg(:,[1,27,28],iSphere)=nan; %scenes 1/28 are cartoons, scene 27 ends the first scan, no event boundary representation in the data
    end
    
    save([pathO files(iSub).name(1:end-4) '_reps.mat'],'repEB','sceneAvg','repP','repF','sphereCenters');
    
end
end

