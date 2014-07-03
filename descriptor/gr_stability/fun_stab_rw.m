function [trkWAll, stab_rw] = fun_stab_rw(sort_stat_trk_ind, curSubdataTrkWAll)

% FUN_STAB_RW: Graph-based group member relative stability
%              Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%%
trkWTimeLine = curSubdataTrkWAll(:,1);

% if there is [] some time, delete it
emptylog = cellfun(@isempty, trkWTimeLine);
trkWTimeLine(emptylog) = [];

trkIndAll = unique(sort_stat_trk_ind(:,2));
trkWAll = zeros(max(trkIndAll), max(trkIndAll), length(trkWTimeLine)); % the whole graph probability matrix from t0 to t1

for t = 1 : length(trkWTimeLine)
    curTrkW = trkWTimeLine{t,1};
    curTrkInd = curTrkW(2:end,1);
    
    trkWAll(curTrkInd, curTrkInd, t) = curTrkW(2:end,2:end);
end

% compute KL-divergence over each point [Tip: if the point disappeared, the row may all zero, in this case, stop]
num = 0;
for i = 1 : size(trkWAll,1)
    curTrkWTimeLine = reshape(trkWAll(i,:,:), size(trkWAll,2), size(trkWAll,3));
    zero_idx = find(sum(abs(curTrkWTimeLine),1)==0, 1);
    if length(find(curTrkWTimeLine==0)) >= (size(curTrkWTimeLine,1)*size(curTrkWTimeLine,2))

    elseif isempty(zero_idx)
        num = num + 1;
        curTrkWnonzero = curTrkWTimeLine;
        for j = 1 : size(curTrkWnonzero,2)-1
            stab(j) = KLDiv([1:size(curTrkWnonzero,1)]',curTrkWnonzero(:,j),curTrkWnonzero(:,j+1),'sym');
        end
        stab_rw(num,1) = mean(stab);
    else
        stab = [];
        for j = 1 : size(curTrkWTimeLine,2)-1
            if sum(abs(curTrkWTimeLine(:,j))) ~= 0 && sum(abs(curTrkWTimeLine(:,j+1))) ~= 0
                stab = [stab KLDiv([1:size(curTrkWTimeLine,1)]',curTrkWTimeLine(:,j),curTrkWTimeLine(:,j+1),'sym')];
            end
        end
        if ~isempty(stab)
            num = num + 1;
            stab_rw(num,1) = mean(stab);
        end
    end
end






