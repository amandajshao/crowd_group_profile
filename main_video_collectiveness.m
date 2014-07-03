%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% main function for collectiveness %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.
clc;clear;close all

%% Descriptor -- Collectiveness (A-fit error within each group)
path = '.\';
path_gr = [path,'result_groupDet_new\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
path_img_dir = xls(2:end,1);

% initialization and parameter setting
group_size_th = 25; 
fit_time_len = 8; % the same as the coherent group detection's parameter

collectiveness(length(path_img_dir),1).coll_gr = []; 
collectiveness(length(path_img_dir),1).coll_gr_mean = []; 
collectiveness(length(path_img_dir),1).gr_len = [];
for file_n = 1 : length(path_img_dir)
    file_name = path_img_dir{file_n};
    fprintf('Group descriptor "Collectiveness" for [%d:%s].\n', file_n, file_name);
    
    %% load collective result from group detection
    load([path_gr, 'trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
    load([path_gr, 'trks_', file_name, '.mat'], 'trks');
    load([path_gr, 'A_1_', file_name, '.mat'], 'A');
    load([path_gr, 'color_1_', file_name, '.mat'], 'color_ind')
    
    %%
    trkClusterNumTime = max(trkClusterTimeLine);
    [trkTime, lenTime, nTrks, trkTimeLine] = fun_trkInfo(trks);
    t_seq = find(trkClusterNumTime ~= 0);
    
    %% Do not need too long time (can be tuned)
    loca = cellfun(@findstr, xls(:,1), repmat({file_name}, size(xls(:,1))), 'UniformOutput', false);
    [t_loc, ~, ~] = find(~cellfun(@isempty, loca) == 1);
    t_start = fun_cell2num(xls(t_loc,5));
    t_end = min(t_seq(end),fun_cell2num(xls(t_loc,6)));

    %%
    cur_trk_ind = find(trkClusterTimeLine(:,t_start)~=0);
    cur_gr_ind = trkClusterTimeLine(cur_trk_ind,t_start);
    data = fun_curX(trks, nTrks, trkTime, t_start, cur_trk_ind);
    [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
    clusterValue = unique(cur_gr_ind);
    
    %%
    rmse_trk_mem = []; group_size = [];
    for grSele = 1 : length(clusterValue)
        clusterV = clusterValue(grSele);
        for curTime = t_start : t_end
            t_oo = find(t_seq==curTime);
            if isempty(t_oo)
                continue;
            end
            A_cur = A(:,t_oo);
            cur_color_ind = color_ind{t_oo};
            cur_color_ind_cluster = find(cur_color_ind==clusterV,1);
            
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
            
            %%
            trk_mem_cur = [];
            if ~isempty(cur_color_ind_cluster)
                group_size{clusterV,curTime} = size(subdata,1);
                A_cur_cluster = A_cur{cur_color_ind_cluster, 1};
                rmse_trk_mem_cur = zeros(size(ind,1),1);
                for i = 1 : size(ind,1)
                    trk_mem = trks(cur_trk_ind(ind(i)));
                    trk_mem_t = find(trk_mem.t >= curTime & trk_mem.t <= curTime+fit_time_len-1);
                    trk_mem_cur = [trk_mem.x(trk_mem_t)'; trk_mem.y(trk_mem_t)'; ones(1, length(trk_mem_t))];
                    
                    % use A_cur_cluster to generate trk_mem_cur_gen
                    trk_mem_cur_gen = ones(3, size(trk_mem_cur,2));
                    trk_mem_cur_gen(:,1) = trk_mem_cur(:,1);
                    for fit_t = 2 : size(trk_mem_cur,2)
                        trk_mem_cur_gen(:,fit_t) = A_cur_cluster * trk_mem_cur_gen(:,fit_t-1);
                    end
                    
                    % compute fit error -- RMSE
                    rmse_trk_mem_cur(i,1) = sqrt(sum(sum((trk_mem_cur_gen-trk_mem_cur).^2)) / size(trk_mem_cur,2));
                end
                
                if isempty(ind)
                    rmse_trk_mem{clusterV,curTime} = [];
                else
                    rmse_trk_mem{clusterV,curTime} = sum(rmse_trk_mem_cur)/length(rmse_trk_mem_cur);
                end
            end
        end
    end

    
    %% record: only groupSize larger than threshold
    group_size_mean = (sum(cellfun(@sum, group_size),2))./(sum(~cellfun(@isempty,group_size),2));
    rmse_trk_mem_mean = (sum(cellfun(@sum, rmse_trk_mem),2))./(sum(~cellfun(@isempty,rmse_trk_mem),2));
    group_len = sum(~cellfun(@isempty,rmse_trk_mem),2);
    
    count1 = []; count2 = []; grSize = [];
    for gr_n = 1 :length(group_size_mean)
        if group_size_mean(gr_n) > group_size_th
            collectiveness(file_n).coll_gr = [collectiveness(file_n).coll_gr; rmse_trk_mem(gr_n,:)];
            collectiveness(file_n).coll_gr_mean = [collectiveness(file_n).coll_gr_mean; rmse_trk_mem_mean(gr_n,:)];
            collectiveness(file_n).gr_len = [collectiveness(file_n).gr_len; group_len(gr_n,:)];
            grSize = [grSize, group_size_mean(gr_n)];
            % for ms
            rmse_temp = cellfun(@sum,rmse_trk_mem(gr_n,:));
            rmse_temp = rmse_temp(rmse_temp~=0);
            temp1 = round(length(rmse_temp)/2);
            count1 = [count1; length(find(rmse_temp(2:temp1)-rmse_temp(1:temp1-1) > 0)) / (temp1-1+eps)];
            temp2 = length(rmse_temp) - temp1;
            count2 = [count2; length(find(rmse_temp(temp1+1:length(rmse_temp))-rmse_temp(temp1:length(rmse_temp)-1) > 0)) / (temp2+eps)];
        end
    end
    
    if length(collectiveness(file_n).coll_gr_mean) == 1
        collectiveness(file_n).coll_v_mean = collectiveness(file_n).coll_gr_mean;
        collectiveness(file_n).coll_v_var = 0;
        
        count = [count1, count2];
        % for v_ms
        [temp1_val, temp1] = max(abs(count(:,1)./count(:,2)));
        [temp2_val, temp2] = max(abs(count(:,2)./count(:,1)));
        if temp1_val > temp2_val
            collectiveness(file_n).coll_count = count(temp1,:);
        elseif temp1_val < temp2_val
            collectiveness(file_n).coll_count = count(temp2,:);
        else
            collectiveness(file_n).coll_count = count(temp1,:);
        end
    elseif length(collectiveness(file_n).coll_gr_mean) > 1
        collectiveness(file_n).coll_v_mean = mean(collectiveness(file_n).coll_gr_mean);
        collectiveness(file_n).coll_v_var = std(collectiveness(file_n).coll_gr_mean);
        
        count = [count1, count2];
        % for v_ms
        [temp1_val, temp1] = max(abs(count(:,1)./count(:,2)));
        [temp2_val, temp2] = max(abs(count(:,2)./count(:,1)));
        if temp1_val > temp2_val
            collectiveness(file_n).coll_count = count(temp1,:);
        elseif temp1_val < temp2_val
            collectiveness(file_n).coll_count = count(temp2,:);
        else
            collectiveness(file_n).coll_count = count(temp1,:);
        end
    end
    collectiveness(file_n).file_name = file_name;
    fprintf('Done!\n');
end

path_result = [path, 'result_groupDescr_new\'];
mkdir(path_result);
save([path_result,'collectiveness.mat'], 'collectiveness')



