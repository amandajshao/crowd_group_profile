function [groupPoint, group_size] = fun_conflict_val(trkClusterTimeLine, A, trks, color_ind, ...
    t_0, t_1, t_seq, nTrks, trkTime, K_percent, clusterValue, grSele, ...
    gen_mem_time_len, fit_time_len, fric_norm)

% FUN_CONFLICT_VAL: Summary of this function goes here
%                   Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

groupPoint = [];

group_size = [];
for curTime = t_0 : t_1
    t_oo = find(t_seq==curTime);
    if isempty(t_oo)
        continue;
    end
    A_cur = A(:,t_oo);
    cur_color_ind = color_ind{t_oo};
    trk_mem_vel = []; fric_mem = []; fric_num_mem = [];
    
    % prepare data
    cur_trk_ind = find(trkClusterTimeLine(:,curTime)~=0);
    cur_gr_ind = trkClusterTimeLine(cur_trk_ind,curTime);
    data = fun_curX(trks, nTrks, trkTime, curTime, cur_trk_ind);
    
    % preprocess data
    [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
    ind = find(cur_gr_ind==clusterValue(grSele));
    trk_mem = cur_trk_ind(ind); 
    if isempty(trk_mem)
        cur_trk_ind = [];
        cur_gr_ind = [];
        data = [];
        break;
    end
    subdata = data(ind,:);
    group_size{1,curTime} = size(subdata,1);
        
    % K need to be corrsponding to group size, and changes as group size changes over time
    K(curTime) = ceil(size(subdata,1) * K_percent);    
    dist_data = distmat(data);
    [~,NNIndex] = gacMink(dist_data, K(curTime)+1, 2);
    
    % compute friction value -- use A of nn to fit member, error is the friction value
    cur_color_ind_cluster = find(cur_color_ind==clusterValue(grSele), 1);
    
    if ~isempty(cur_color_ind_cluster) % this if is for disappeared clusterValue(clusterSelected)
        A_cur_cluster = A_cur{cur_color_ind_cluster, 1};
        for i = 1 : size(subdata,1)
            trk_mem_gen = []; rmse_trk_nn_cur_gen = [];
            knn_inds = NNIndex(ind(i), :)';
            trk_mem_cur = trks(cur_trk_ind(knn_inds(1)));
            knn_inds = knn_inds(2:end); % knn_inds(1) is subdata(i)
            knn_cluster = cur_gr_ind(knn_inds);
            i_temp = 0;
            % compute each member's velocity in order to calibration later
            trk_mem_gen(1,1) = trk_mem_cur.x(trk_mem_cur.t == curTime);
            trk_mem_gen(2,1) = trk_mem_cur.y(trk_mem_cur.t == curTime);
            trk_mem_gen(3,1) = 1;
            for j = 2 : gen_mem_time_len
                trk_mem_gen(:,j) = A_cur_cluster * trk_mem_gen(:,j-1);
            end
            trk_mem_cur_vel = (trk_mem_gen(:,end)-trk_mem_gen(:,1))/(gen_mem_time_len-1);
            trk_mem_vel(i,:) = trk_mem_cur_vel(1:2)';
            for j = 1 : length(knn_inds)
                trk_nn_cur = []; 
                if knn_cluster(j) ~= clusterValue(grSele) && ~isempty(intersect(knn_cluster(j), clusterValue))
                    % this knn is from a different cluster, use A_cur_cluster to fit it
                    i_temp = i_temp + 1;
                    trk_nn = trks(cur_trk_ind(knn_inds(j)));
                    trk_nn_t = find(trk_nn.t >= curTime & trk_nn.t <= curTime+fit_time_len-1);
                    trk_nn_cur(1,:) = trk_nn.x(trk_nn_t);
                    trk_nn_cur(2,:) = trk_nn.y(trk_nn_t);
                    trk_nn_cur(3,:) = ones(1, length(trk_nn_t));
                    % use A_cur_cluster to generate trk_nn_cur_gen
                    trk_nn_cur_gen = ones(3, size(trk_nn_cur,2));
                    trk_nn_cur_gen(:,1) = trk_nn_cur(:,1);
                    for fit_t = 2 : size(trk_nn_cur,2)
                        trk_nn_cur_gen(:,fit_t) = A_cur_cluster * trk_nn_cur_gen(:,fit_t-1);
                    end
                    % compute fit error -- RMSE
                    rmse_trk_nn_cur_gen(1,i_temp) = sqrt(sum(sum((trk_nn_cur_gen-trk_nn_cur).^2))) / size(trk_nn_cur,2);
                end
            end
            if i_temp == 0
                % that is, this subdata(i) does not have any friction
                fric_mem(i) = 0;
                fric_num_mem(i) = i_temp;
            else
                if fric_norm == 1
                    fric_mem(i) = sum(rmse_trk_nn_cur_gen)/i_temp;
                else
                    fric_mem(i) = sum(rmse_trk_nn_cur_gen);
                end
                fric_num_mem(i) = i_temp;
            end
        end
        % store friction info : [curTime x y v_x v_y fricValue fricNum]
        groupPoint = [groupPoint; [ones(size(subdata,1),1)*curTime, subdata, trk_mem_vel, fric_mem', fric_num_mem']];
    end
end


