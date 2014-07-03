function [CF_init, group_ref, id_ref, color_ind, A, t, trkClusterTimeLine, trks] = ...
    cohGr_det(trks, path_img, frames, paraSet)

% COHGR_DET: Core function for coherent group detection
%            Detailed explanation goes here
% --------------------------------------------------------------------- %
% Input: trks            -- klt tracking results
%        path_img        -- image sequence folder
%        paraSet         -- parameter set
%        file_name       -- image sequence name
%        lenTime_def     -- the length of frames to detect group (predefined or not)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate/Output variable:
%       CF_init      -- the indices of coherent trks at every frame
%       group_init   -- the initial group trks within t~t+t_interval-1 from CF itself
%       group_seed   -- the reliable seeddate group trks within t~t+t_interval-1
%       group_ref    -- the A-fit-based group refinement
%       color_ind    -- color index association for group visualization
%       A            -- collective transition prior
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% "num" is help to determine the ones belong to one group
for i = 1 : length(trks)
    trks(i).num = 0;
    trks(i).num_temp = 0;
    % mark the trk, if the track has been completed or been used, mark=1
    trks(i).mark = 0; 
end

%% Initialization
len_fr_def = paraSet.len_fr;
kf_max_iter = paraSet.max_iter;
kf_t_intv = paraSet.t_intv_kf;
gr_err_th = paraSet.err_th;
gr_t_intv = paraSet.t_intv;
gr_dv_th_seed = paraSet.dv_th_seed;
color_all = [1:100]';
[trkTime, lenTime, nTrks, ~] = fun_trkInfo(trks);
lenTime = min(lenTime,length(frames));

%% Step1 : Coherent filtering over the video, and store all the group index in "CF_init"
if len_fr_def ~= 0 % if the length of frames is predefined
    fprintf('Coherent filtering begins...\n');
    [CF_init,trks] = fun_CF(trks, path_img, 1, min(len_fr_def,lenTime), paraSet);
else
    fprintf('Coherent filtering begins...\n');
    [CF_init,trks] = fun_CF(trks, path_img, 1, lenTime, paraSet);
end

%% Step2-4 : Update group every t_interval time
trkClusterTimeLine(nTrks,max(trkTime(:))) = uint8(0);
t_all = 0;
t_cur = paraSet.start_fr;
t_stop = length(CF_init)-gr_t_intv+1; % can set as other number to control detect groups within t_stop frames
fprintf('CT-based group detection begins...\n');
while t_cur <= t_stop
    fprintf('Doing frame %d.\n', t_cur);
    if isempty(CF_init{1,t_cur}) % when CF break sometime, let t=t+1
%         fun_grVisual([], path_img, frames, t_cur); % visualization
        t_all = t_all + 1;
        group_ref{1,t_all} = []; id_ref{1,t_all} = []; color_ind{1,t_all} = []; A{1,t_all} = []; t{1,t_all} = 0;
        t_cur = t_cur + 1;
    else        
        len_CF = zeros(length(CF_init{1,t_cur}.gr),1);
        for i = 1 : length(CF_init{1,t_cur}.gr)
            len_CF(i,1) = len_CF(i,1) + length(CF_init{1,t_cur}.gr{i,1}(:,1));
        end
        
        gr_num = 0; % group number for our detection
        g = 1; % group number for CF_init
        t_all = t_all + 1; % the t_all interval
        while sum(len_CF)>3*length(len_CF) && g<=length(len_CF) % add this while to deal with those not in any group
            gr_num = gr_num + 1;
            
            %% Step2 : identify anchor tracklet
            [trkAnchor_ind_cur, trkAnchor_cur] = fun_anchorTrk(CF_init{1,t_cur}.gr{g,:}, trks);
            % record anchor tracklets and their indices
            trkAnchor_ind{gr_num,t_all} = trkAnchor_ind_cur; trkAnchor{gr_num,t_all} = trkAnchor_cur;
            
            %% Step3 : discover seeding tracklets
            if isempty(trkAnchor_ind_cur) % do not find anchor tracklet, use another condition to select tracklet
                trkAnchor{gr_num,t_all} = [];
                trkAnchor_ind{gr_num,t_all} = [];
                gr_num = gr_num - 1;
                g = g + 1;
                continue;
            else
                %% Step3_1 : detect initial group (group_init) (this is from pure CF itself without any process)
                %%           and discover seeding group (group_seed)
                [group_init_cur, group_seed_cur, id_init_cur, id_seed_cur] = ...
                    fun_groupInit(CF_init, trks, trkAnchor_ind_cur, trkAnchor_cur, t_cur, gr_t_intv, gr_dv_th_seed);
                % record intial group and seeding group
                group_init{gr_num,t_all} = group_init_cur; group_seed{gr_num,t_all} = group_seed_cur; id_init{gr_num,t_all} = id_init_cur; id_seed{gr_num,t_all} = id_seed_cur;
                
                if length(group_seed{gr_num,t_all}) < 3 % few tracklets in seeding group, redetection
                    %% Reselect the "trkTrack"
                    CF_init{1,t_cur}.gr{g,1}(CF_init{1,t_cur}.gr{g,1}(:,1)==trkAnchor_ind_cur,2) = -1;
                    group_init{gr_num,t_all} = [];
                    group_seed{gr_num,t_all} = [];
                    id_init{gr_num,t_all} = [];
                    id_seed{gr_num,t_all} = [];
                    trkAnchor{gr_num,t_all} = [];
                    trkAnchor_ind{gr_num,t_all} = [];
                    gr_num = gr_num - 1;
                else
                    %% Step3_2 : Kalman filter training on seeding group
                    A_cur = fun_kf(group_seed_cur, kf_max_iter, t_cur, kf_t_intv);
                    A{gr_num,t_all} = A_cur;
                    
                    %% Step4 : A-fit-based refine group
                    [group_ref_cur, id_ref_cur, trkLatent_cur] = fun_groupRef(trkAnchor_cur, ...
                        group_init_cur, id_init_cur, id_seed_cur, A_cur, t_cur, t_cur+kf_t_intv-1, gr_err_th);
                    % record refined group and latent tracklet
                    group_ref{gr_num,t_all} = group_ref_cur;
                    id_ref{gr_num,t_all} = id_ref_cur;
                    trkLatent{gr_num,t_all} = trkLatent_cur;                    
                    % add a condition to judge that we add this group or not
                    if length(group_ref{gr_num,t_all}) < 3 % do not add too small group
                        CF_init{1,t_cur}.gr{g,1}(CF_init{1,t_cur}.gr{g,1}(:,1)==trkAnchor_ind_cur,2) = -1;
                        group_init{gr_num,t_all} = [];
                        group_seed{gr_num,t_all} = [];
                        group_ref{gr_num,t_all} = [];
                        id_init{gr_num,t_all} = [];
                        id_seed{gr_num,t_all} = [];
                        id_ref{gr_num,t_all} = [];
                        trkAnchor{gr_num,t_all} = [];
                        trkAnchor_ind{gr_num,t_all} = [];
                        A{gr_num,t_all} = [];
                        gr_num = gr_num - 1;
                    end
                end
            end
            
            %% Delete trk in "group_ref" from "CF_init" (their indices are all set to zero)
            if gr_num ~= 0
                [CF_init, len_CF] = fun_CFUpd(CF_init, t_cur, id_ref_cur);
            else
                [~, len_CF] = fun_CFUpd(CF_init, t_cur);
            end
            if len_CF(g,1) <= 3
                g = g + 1;
            end
        end
                
            
        %% visualize: the groups within this time interval by time
        if gr_num~= 0
            %% Step5 : Update color_ind for each group in this time interval
            if t_all == 1 % the first time interval
                color_ind{:,t_all} = [1:gr_num]';
                color_all = setdiff(color_all, color_ind{:,t_all});
            elseif sum(~cellfun(@isempty, id_ref(:,t_all-1))) == 0
                color_ind{:,t_all} = [1:gr_num]';
                color_all = setdiff(color_all, color_ind{:,t_all});
            else
                [color_ind, color_all] = fun_colorAss(color_ind, color_all, t_all, id_ref, trkLatent);
            end
            % record trkClusterTimeLine
            for i = 1 : sum(~cellfun(@isempty, id_ref(:,t_all)))
                trkClusterTimeLine(id_ref{i,t_all}, t_cur) = uint8(color_ind{1,t_all}(i,1));
            end
%             fun_grVisual(group_ref, path_img, frames, t_cur, t_cur, color_ind, t_all) % visualization
            t{1,t_all} = t_cur;
            t_cur = t_cur + 1;
        else
%             fun_grVisual([], path_img, frames, t_cur) % visualization
            t{1,t_all} = t_cur;
            t_cur = t_cur + 1;
        end
    end
    
end









