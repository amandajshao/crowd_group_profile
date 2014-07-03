function [pairSet,correSet] = CF_neighbor2pair(neighborSet,correlationSet,d)

zeroNeighborSet=neighborSet{1,1};
nPoint=size(zeroNeighborSet,1);


pairSet=[];
correSet=[];

for i=1:nPoint  
    
    curIntersect=zeroNeighborSet(i,:);
    curCorre=zeros(1,size(zeroNeighborSet,2));
    for j=1:d
        nextNeighborSet=neighborSet{1,j};
        nextCorreSet=correlationSet{1,j};
        [curIntersect,ia,ib]=intersect(curIntersect,nextNeighborSet(i,:)); % the intersection of neighborhood
        curCorre=curCorre(1,ia)+nextCorreSet(i,ib);
    end   
   
    if ~isempty(curIntersect)

        
        correSet=[correSet curCorre./d];
        pairSet=[pairSet [i*ones(1,length(curIntersect));curIntersect]];  
    end
end


end

