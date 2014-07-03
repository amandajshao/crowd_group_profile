
function clusterIndex=pair2cluster(pairwiseData,totalNum)

clusterNum=0;
clusterIndex=zeros(1,totalNum);

for i=1:size(pairwiseData,2)
    curPair=pairwiseData(:,i);
    curPairAlabel=clusterIndex(1,curPair(1));
    curPairBlabel=clusterIndex(1,curPair(2));
    
    if curPairAlabel==0 && curPairBlabel==0
        clusterNum=clusterNum+1;
        curPairLabel=clusterNum;
        clusterIndex(1,curPair(1))=curPairLabel;
        clusterIndex(1,curPair(2))=curPairLabel;
    elseif curPairAlabel~=0 && curPairBlabel==0
        clusterIndex(1,curPair(2))=curPairAlabel;
    elseif curPairBlabel~=0 && curPairAlabel==0
        clusterIndex(1,curPair(1))=curPairBlabel;
    else
        combineLabel=min(curPairAlabel,curPairBlabel);
        clusterIndex(1,find(clusterIndex==curPairAlabel))=combineLabel;
        clusterIndex(1,find(clusterIndex==curPairBlabel))=combineLabel;
    end
    

end

newClusterNum=0;
for i=1:max(clusterNum)
    curClusterIndex=find(clusterIndex==i);
    if length(curClusterIndex)<5
        clusterIndex(curClusterIndex)=0;
   
    else
        newClusterNum=newClusterNum+1;
        clusterIndex(curClusterIndex)=newClusterNum;
    end
end
    

end
