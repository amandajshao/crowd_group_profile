function [group_ref_cur, id_ref_cur, trkLatent] = ...
    fun_groupRef(trkAnchor_cur, group_init_cur, id_init_cur, id_seed_cur, A_cur, t_begin, t_end, err_th)

% FUN_GROUPREF: Group refinement
%               Detailed explanation goes here
% --------------------------------------------------------------------- %
% Input :  CF_init       -- the indices of coherent trks at every frame
%          trkAnchor     -- anchor tracklet (use this to detect seeding group)
% Output : group_ref_cur -- refined group with small fit error
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% check each trk based on A-fit method
num_temp = 0;
for i = 1 : length(group_init_cur)
    group_gt = [];
    group_gen = [];
    
    % use A fit each CF_trk within the period between time_begin and time_end
    group_init_cur(i).x = group_init_cur(i).x(group_init_cur(i).t>=t_begin & group_init_cur(i).t<=t_end);
    group_init_cur(i).y = group_init_cur(i).y(group_init_cur(i).t>=t_begin & group_init_cur(i).t<=t_end);
    group_init_cur(i).t = group_init_cur(i).t(group_init_cur(i).t>=t_begin & group_init_cur(i).t<=t_end);
    group_gen(:,1) = [group_init_cur(i).x(1); group_init_cur(i).y(1); 1];
    group_gt(:,1) = [group_init_cur(i).x(1); group_init_cur(i).y(1); 1];
    for j = 2 : length(group_init_cur(i).t)
        group_gen(:,j) = A_cur * group_gen(:,j-1);
        group_gt(:,j) = [group_init_cur(i).x(j); group_init_cur(i).y(j); 1];
    end
    
    % compute fit error -- RMSE
    mse_CF_trk = sqrt(sum(sum((group_gen-group_gt).^2))) / length(group_init_cur(i).t);
%     fprintf('mse_CF_trk:%f\n', mse_CF_trk);
    if mse_CF_trk < err_th
        num_temp = num_temp + 1;
        group_ref_cur(num_temp,1).x = group_init_cur(i).x;
        group_ref_cur(num_temp,1).y = group_init_cur(i).y;
        group_ref_cur(num_temp,1).t = group_init_cur(i).t;
        id_ref_cur(num_temp,1) = id_init_cur(i,1);
    elseif ~isempty(intersect(id_init_cur(i,1),id_seed_cur))
        num_temp = num_temp + 1;
        group_ref_cur(num_temp,1).x = group_init_cur(i).x;
        group_ref_cur(num_temp,1).y = group_init_cur(i).y;
        group_ref_cur(num_temp,1).t = group_init_cur(i).t;
        id_ref_cur(num_temp,1) = id_init_cur(i,1);
    end
end

%% 
t_temp = find(trkAnchor_cur.t>=t_begin & trkAnchor_cur.t<=t_end);
trkLatent(1,:) = trkAnchor_cur.x(t_temp)';
trkLatent(2,:) = trkAnchor_cur.y(t_temp)';
trkLatent(3,:) = ones(1,length(t_temp));
trkLatent(4,:) = trkAnchor_cur.t(t_temp)';


end