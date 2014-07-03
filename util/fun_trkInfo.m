function [trkTime, lenTime, nTrks, trkTimeLine] = fun_trkInfo(trks)

% FUN_TRKINFO: Summary of this function goes here
%              Detailed explanation goes here
% --------------------------------------------------------------------- %
% Output: trkTime -- store all trks' start_frame and end_frame
%         lenTime -- max frame
%         nTrks   -- the number of trks
%         trkTimeLine -- record each track's period
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.


nTrks = length(trks);
trkTime = zeros(2,nTrks);
for i = 1 : length(trks)
    trkTime(:,i) = [trks(i).t(1); trks(i).t(end)];
end
lenTime = max(trkTime(2,:));

trkTimeLine = zeros(nTrks,max(trkTime(:)));
for i = 1 : nTrks
    trkTimeLine(i,trkTime(1,i):trkTime(2,i)) = 1;
end

end