function [color_ind, color_all] = fun_colorAss(color_ind, color_all, t_all, id_ref, trkLatent)

% FUN_COLORASS: Associate color_ind for each group in this time interval
%               Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% in the following intervals, compare with the previous 2 frames to determine the color for each group
for i = 1 : sum(~cellfun(@isempty, id_ref(:,t_all)))
    common_id = [];
    for j = 1: sum(~cellfun(@isempty, id_ref(:,t_all-1)))
        common_id(j,1) = length(intersect(id_ref{i,t_all}, id_ref{j,t_all-1}));
    end
    [v_temp, id_temp] = max(common_id(:,1));
    if v_temp == 0
        if t_all == 2
            if i == 1
                color_ind{:,t_all}(i,1) = min(setdiff(color_all, color_ind{:,t_all-1}));
            else
                color_ind{:,t_all}(i,1) = min(setdiff(setdiff(color_all, color_ind{:,t_all-1}), color_ind{:,t_all}));
            end
        else
            for j = 1: sum(~cellfun(@isempty, id_ref(:,t_all-2)))
                common_id(j,1) = length(intersect(id_ref{i,t_all}, ...
                    id_ref{j,t_all-2}));
            end
            [v_temp, id_temp] = max(common_id(:,1));
            if v_temp == 0
                if i == 1
                    color_ind{:,t_all}(i,1) = min(setdiff(color_all, color_ind{:,t_all-2}));
                else
                    color_ind{:,t_all}(i,1) = min(setdiff(setdiff(color_all, color_ind{:,t_all-2}), color_ind{:,t_all}));
                end
            else
                x_temp1=trkLatent{i,t_all}(1,:);
                y_temp1=trkLatent{i,t_all}(2,:);
                x_temp2=trkLatent{id_temp,t_all-2}(1,:);
                y_temp2=trkLatent{id_temp,t_all-2}(2,:);
                
                direct_temp = dot([x_temp1(end)-x_temp1(1),y_temp1(end)-y_temp1(1)], ...
                    [x_temp2(end)-x_temp2(1),y_temp2(end)-y_temp2(1)]) / ...
                    (norm([x_temp1(end)-x_temp1(1),y_temp1(end)-y_temp1(1)]) * ...
                    norm([x_temp2(end)-x_temp2(1),y_temp2(end)-y_temp2(1)]));
                
                if direct_temp > 0.4
                    color_ind{:,t_all}(i,1) = color_ind{:,t_all-2}(id_temp,1);
                else
                    if i == 1
                        color_ind{:,t_all}(i,1) = min(setdiff(color_all, color_ind{:,t_all-2}));
                    else
                        color_ind{:,t_all}(i,1) = min(setdiff(setdiff(color_all, color_ind{:,t_all-2}), color_ind{:,t_all}));
                    end
                end
            end
        end
    else
        x_temp1=trkLatent{i,t_all}(1,:);
        y_temp1=trkLatent{i,t_all}(2,:);
        x_temp2=trkLatent{id_temp,t_all-1}(1,:);
        y_temp2=trkLatent{id_temp,t_all-1}(2,:);
        direct_temp = dot([x_temp1(end)-x_temp1(1),y_temp1(end)-y_temp1(1)], ...
            [x_temp2(end)-x_temp2(1),y_temp2(end)-y_temp2(1)]) / ...
            (norm([x_temp1(end)-x_temp1(1),y_temp1(end)-y_temp1(1)]) * ...
            norm([x_temp2(end)-x_temp2(1),y_temp2(end)-y_temp2(1)]));
        if direct_temp > 0.4
            color_ind{:,t_all}(i,1) = color_ind{:,t_all-1}(id_temp,1);
        else
            if i == 1
                color_ind{:,t_all}(i,1) = min(setdiff(color_all, color_ind{:,t_all-1}));
            else
                color_ind{:,t_all}(i,1) = min(setdiff(setdiff(color_all, color_ind{:,t_all-1}), color_ind{:,t_all}));
            end
        end
    end
end


end