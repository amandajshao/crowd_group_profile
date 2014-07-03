function [neighborSet,correlationSet] = CF_neighbor(allXset,allVset,K)
%SUP_NEIGHBOR Summary of this function goes here
%   Detailed explanation goes here
d=length(allXset);
nPoint=size(allXset{1,1},2);
neighborSet=cell(1,d);
correlationSet=cell(1,d);

for j=1:d
%      curAllX=allXset{1,j};
     curAllV=allVset{1,j};
     curNeighborGraph=zeros(nPoint,K);
     curCorrelationGraph=zeros(nPoint,K);
   

     for i=1:nPoint
        curAllX=allXset{1,j}(1:2,:);
         
        curX=curAllX(:,i);
        distance=repmat(curX,[1 nPoint])-curAllX;
        distance=sqrt(sum(distance.^2));
        [B,IX] = sort(distance,'ascend');   
        for kNN=1:min(K,numel(IX))-1
        
            if sqrt(sum(curAllV(:,i).^2))>0 && sqrt(sum(curAllV(:,IX(kNN+1)).^2)) %&&B(j+1)<dis_th
            
                coefficient=curAllV(:,i)'*curAllV(:,IX(kNN+1));
                coefficient=coefficient/(sqrt(sum(curAllV(:,i).^2))*sqrt(sum(curAllV(:,IX(kNN+1)).^2)));          
                curNeighborGraph(i,kNN)=IX(kNN+1);
                curCorrelationGraph(i,kNN)=coefficient;
            end       
        end
        
     end
     neighborSet{1,j}=curNeighborGraph;
     correlationSet{1,j}=curCorrelationGraph;

end


