function [group_init_cur, group_seed_cur, id_init_cur, id_seed_cur] = ...
    fun_groupInit(CF_init, trks, trkAnchor_ind, trkAnchor, t, t_interval, dv_th)

% FUN_GROUPINIT: Identify intial group and seeding group
%                Detailed explanation goes here
% --------------------------------------------------------------------- %
% Input :  CF_init    -- the indices of coherent trks at every frame
%          trkAnchor  -- anchor tracklet (use this to detect seeding group)
% Output : group_init -- initial group trks within t~t+t_interval-1 from CF itself
%          group_seed -- seeding group within t~t+t_interval-1
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% find max_trk_ind and min_trk_ind in this time interval
trk_maxInd_temp = []; trk_minInd_temp = [];
for tt = t : t+t_interval-1
    if ~isempty(CF_init{1,tt})
        trk_maxInd_temp = [trk_maxInd_temp; max(cellfun(@(x) max(x(:)), CF_init{1,tt}.gr))];
        trk_minInd_temp = [trk_minInd_temp; min(cellfun(@(x) min(x(:)), CF_init{1,tt}.gr))];
    end
end
trk_maxInd = max(trk_maxInd_temp);
trk_minInd = max(min(trk_minInd_temp),1);

%% group_seed -- find a reliable and robust seeding group
for i = 1 : length(trks)
    trks(i).num_temp = 0;
end
trks(trkAnchor_ind).num_temp = 1000;
num_temp = 0;
for tt = t : t+t_interval-1
    if ~isempty(CF_init{1,tt})
        for g = 1 : length(CF_init{1,tt}.gr(:,1))
            if ~isempty(intersect(CF_init{1,tt}.gr{g,1}(:,1), trkAnchor_ind))
                num_temp = num_temp + 1;
                id_temp = find(CF_init{1,tt}.gr{g,1}(:,1)~=0);
                for i = 1 : length(CF_init{1,tt}.gr{g,1}(id_temp,1))
                    trks(CF_init{1,tt}.gr{g,1}(id_temp(i),1)).num_temp = ...
                        trks(CF_init{1,tt}.gr{g,1}(id_temp(i),1)).num_temp + 1;
                end
                break;
            end
        end
    end
end
[group_seed_cur, id_seed_cur] = gener_group_seed(trks, trkAnchor, trk_minInd, trk_maxInd, ...
    num_temp, t, t_interval, dv_th);

%% group_init -- from pure CF_init itself
group_init_cur = trks(trk_minInd:trk_maxInd)';
id_init_cur = [trk_minInd:trk_maxInd]';
i_temp = 1;
while i_temp <= length(group_init_cur)
    if group_init_cur(i_temp).num_temp == 0
        id_init_cur(i_temp) = [];
        group_init_cur(i_temp) = [];
        i_temp = i_temp - 1;
    end
    i_temp = i_temp + 1;
end

end

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

function [group_seed_cur, id_seed_cur] = gener_group_seed(trks, trkAnchor, trk_minInd, ...
    trk_maxInd, num, t, t_interval, dv_th)

group_seed_cur = trks(trk_minInd:trk_maxInd)';
id_seed_cur = [trk_minInd:trk_maxInd]';
trkInterval = find(trkAnchor.t>=t & trkAnchor.t<=t+t_interval-1);
trkSec.x = trkAnchor.x(trkInterval);
trkSec.y = trkAnchor.y(trkInterval);
trkSec.t = trkAnchor.t(trkInterval);

%% bigger "num_temp" means that this coherent track is stable
i = 1;
while i <= length(group_seed_cur)
    groupInterval = find(group_seed_cur(i).t>=t & group_seed_cur(i).t<=t+t_interval-1);
    groupSec.x = group_seed_cur(i).x(groupInterval);
    groupSec.y = group_seed_cur(i).y(groupInterval);
    groupSec.t = group_seed_cur(i).t(groupInterval);
    groupSec.num_temp = group_seed_cur(i).num_temp;
    
    [v,i1,i2] = intersect(trkSec.t, groupSec.t);
    if isempty(v)
        id_seed_cur(i) = [];
        group_seed_cur(i) = [];
        i = i - 1;
    else    
        temp1 = [trkSec.x(i1(end))-trkSec.x(i1(1)), trkSec.y(i1(end))-trkSec.y(i1(1))];
        temp2 = [groupSec.x(i2(1))-trkSec.x(i1(1)), groupSec.y(i2(1))-trkSec.y(i1(1))];
        dv = abs(dot(temp1,temp2) / (norm(temp1)));
        if groupSec.num_temp < ceil(num/3*1) % short time to be together with anchor tracklet
            id_seed_cur(i) = [];
            group_seed_cur(i) = [];
            i = i - 1;
        else % use velocity (direction) rule to determine more accurate and stable coherent track
            if dot([trkSec.x(i1(end))-trkSec.x(i1(1)), trkSec.y(i1(end))-trkSec.y(i1(1))], ...
                    [groupSec.x(i2(end))-groupSec.x(i2(1)), groupSec.y(i2(end))-groupSec.y(i2(1))])/...
                    (norm([trkSec.x(i1(end))-trkSec.x(i1(1)), trkSec.y(i1(end))-trkSec.y(i1(1))]) * ...
                    norm([groupSec.x(i2(end))-groupSec.x(i2(1)), groupSec.y(i2(end))-groupSec.y(i2(1))]))...
                    < dv_th
                id_seed_cur(i) = [];
                group_seed_cur(i) = [];
                i = i - 1;
            end
        end
    end
    i = i + 1;
end

end

