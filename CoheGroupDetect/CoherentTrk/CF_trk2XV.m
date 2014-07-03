function [allXset,allVset] = CF_trk2XV(trks,curTime,d)
%SUP_TRK2XV Summary of this function goes here
%   Detailed explanation goes here
% sampleTrk=trks(1,1);
% nDimension=length(fieldnames(sampleTrk))-1;

allXset=cell(1,d);
allVset=cell(1,d);

for i=1:d
    allXset{1,i}=[];
    allVset{1,i}=[];
end

for i=1:length(trks)
    curStart=trks(i).t(1);curEnd=trks(i).t(end);
    if curTime>curStart && curEnd>=(curTime+d) && length(trks(i).x)>=(curTime+d-curStart)
%         if nDimension==2
            
            curX=[trks(i).x((curTime-curStart):(curTime-curStart)+d-1)';trks(i).y((curTime-curStart):(curTime-curStart)+d-1)'];
            curV=[trks(i).x((curTime-curStart)+1:(curTime-curStart)+d)' - trks(i).x((curTime-curStart):(curTime-curStart)+d-1)';...
                trks(i).y((curTime-curStart)+1:(curTime-curStart)+d)' - trks(i).y((curTime-curStart):(curTime-curStart)+d-1)'];
%         else
%             
%             curX=[trks(i).x((curTime-curStart):(curTime-curStart)+d-1)';trks(i).y((curTime-curStart):(curTime-curStart)+d-1)';trks(i).z((curTime-curStart):(curTime-curStart)+d-1)'];
%             curV=[trks(i).x((curTime-curStart)+1:(curTime-curStart)+d)'-trks(i).x((curTime-curStart):(curTime-curStart)+d-1)';trks(i).y((curTime-curStart)+1:(curTime-curStart)+d)'-trks(i).y((curTime-curStart):(curTime-curStart)+d-1)';trks(i).z((curTime-curStart)+1:(curTime-curStart)+d)'-trks(i).z((curTime-curStart):(curTime-curStart)+d-1)'];        
%         end
        for j=1:d           
            curX_temp = [curX(:,j); i];
            allXset{1,j}=[allXset{1,j},curX_temp];
            allVset{1,j}=[allVset{1,j},curV(:,j)];
        end
        
    end
end


end

