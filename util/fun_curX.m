function [data] = fun_curX(trks, nTrks, trkTime, curTime, cur_trkInd)
% FUN_CURX: Summary of this function goes here
%           Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

curXset = zeros(2,nTrks);
curIndex = find(trkTime(1,:)<=curTime & trkTime(2,:)>=curTime)';
for i = 1 : length(curIndex)
    curTrk = trks(1,curIndex(i));
    pointIndex = find(curTrk.t == curTime);
    curXset(1,curIndex(i)) = curTrk.x(pointIndex);
    curXset(2,curIndex(i)) = curTrk.y(pointIndex);
end
data = curXset(:,cur_trkInd)';

end

