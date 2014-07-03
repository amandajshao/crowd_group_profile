function [trkAnchor_ind, trkAnchor] = fun_anchorTrk(gr_temp, trks)

% FUN_ANCHORTRK: Identify anchor tracklet
%                Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

gr_ind_temp = gr_temp(gr_temp(:,2)~=-1,1); % if trk == -1, that means it cannot be anchor trk
len_trk_temp = [];
for i = 1 : length(gr_ind_temp)
    len_trk_temp(i) = length(trks(gr_ind_temp(i)).x);
end
[~,max_ind] = max(len_trk_temp);
trkAnchor_ind = gr_ind_temp(max_ind,1);   % anchor tracklet index
trkAnchor = trks(trkAnchor_ind); % anchor tracklet

end
