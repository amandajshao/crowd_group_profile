function [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind)
% FUN_CURX_PREPROCESS: Summary of this function goes here
%                      Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% preprocess 1
cur_gr_ind(data(:,1)==0) = [];
cur_trk_ind(data(:,1)==0) = [];
data(data(:,1)==0,:) = [];

% preprocess 2: delete duplicate data
[data_temp,data_ind_1,data_ind_2] = unique(data,'rows','stable');
if length(data_ind_1) ~= length(data_ind_2)
    data = data_temp;
    cur_gr_ind = cur_gr_ind(data_ind_1);
    cur_trk_ind = cur_trk_ind(data_ind_1);
end

end
