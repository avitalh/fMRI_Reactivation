% create BOW similarity matrix
clear
% load lowered-case text, cleaned from punctuations, with empty spaces denoting the two cartoons at scenes 1/28
str = extractFileText('cleanFullText.txt');
textData = split(str,newline);
textData(end)=[]; % last row contains ""

documents = tokenizedDocument(textData);
bag = bagOfWords(documents);

% remove stop words (note that the list of stop words changes between Matlab versions)
newBag = removeWords(bag,stopWords);

% have a look
tbl = topkwords(newBag,20)

% get word counts
data = full(newBag.Counts);

%transform to probabilities
data=(data'./nansum(data'))'; 

% claculate scene distances
for i=1:size(data,1)
    for j=1:size(data,1)
      [JSD(i,j)]=invJSD(data(i,:),data(j,:));
    end
end
imagesc(JSD)

% create a matrix that holds the minimal number of words between each 2 scenes
data(data~=0)=1;
numWords=sum(data');
for i=1:size(data,1)
    for j=1:size(data,1)
      [minWord(i,j)]=min(numWords(i),numWords(j));
    end
end

% ignore cartoon scenes
minWord([1,28],:)=nan;
minWord(:,[1,28])=nan;
