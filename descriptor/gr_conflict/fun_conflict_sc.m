function conflict = fun_conflict_sc(groupPoint, cali_para, sc_para, gv_flag)

% FUN_CONFLICT_SC: Summary of this function goes here
%                  Detailed explanation goes here
% --------------------------------------------------------------------- %
% Input: groupPoint -- [curTime x y conf_value conf_num]
%        gr_n       -- current group index number
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.


%% store all conftion cases
conflict.conflict = groupPoint;

%% initialization
[~,temp1,~] = unique(groupPoint(:,2:3),'rows');
groupPoint = groupPoint(temp1,:);

%% calibration
smooth_size = cali_para.smooth_size;
smooth_sigma = cali_para.smooth_sigma;
pattern_scale = cali_para.pattern_scale;
contour_th = cali_para.contour_th;
[~, points_trans_normalize, density_pattern] = fun_calibration(groupPoint, ...
    smooth_size, smooth_sigma, pattern_scale, gv_flag);
pattern_contour = edge(density_pattern > contour_th);    % get the contour
conflict.density = density_pattern;

%% shape context parameters
nbins_theta = sc_para.nbins_theta;
nbins_r = sc_para.nbins_r;
r_inner = sc_para.r_inner;
r_outer = sc_para.r_outer;
mean_dist_global = sc_para.mean_dist_global;

%% contour
[pc_X, pc_Y] = ind2sub([121,121], find(pattern_contour>0));
Xk1 = [pc_Y, pc_X]; % contour points
contour_val = ones(size(Xk1,1),1);

%% conflict
conf_select_th = ceil(mean(groupPoint(groupPoint(:,6)'>0,6)));
conf_trans = points_trans_normalize(:,groupPoint(:,6)'>=conf_select_th);
conf_val = groupPoint(groupPoint(:,6)>=conf_select_th,6);
Xk2 = conf_trans'; % conflict points

%% shape context over conftional+contour points (conftion points' distribution)
if ~isempty(Xk2) % Xk2:conf; Xk1:contour
    Xk3 = [Xk1; (Xk2.*pattern_scale) + smooth_size + round(pattern_scale/2)];
    Xk3_mark = [zeros(size(Xk1,1),1); ones(size(Xk2,1),1)];
    nsample3 = length(Xk3(:,1));
    out_vec3 = zeros(1, nsample3);
    conf_mark = ones(size(Xk3,1),size(Xk3,1));
    conf_mark(:,end-size(Xk2,1)+1:end) = 0;
    
    for conf_n = 1 : size(Xk2,1)
        Xk3_each = [Xk1; Xk2(conf_n,:).*pattern_scale + smooth_size + round(pattern_scale/2)];
        nsample_each = length(Xk3_each(:,1));
        out_vec_each = zeros(1, nsample_each);
        
        [BH_each,~] = fun_sc_compute(Xk3_each',zeros(1,nsample_each),mean_dist_global,...
            nbins_theta,nbins_r,r_inner,r_outer,out_vec_each);
        BH3_conf(conf_n,:) = BH_each(end,:);
    end
        
    % translate BH3 into histogram form
    BH3_hist = sum(BH3_conf);
    BH3_hist = BH3_hist/size(BH3_conf,2);
    BH3_hist_norm = BH3_hist./sum(BH3_hist);
    
    conflict.sc_conf_contour = BH3_conf;
    conflict.sc_conf_contour_hist = BH3_hist;
    conflict.sc_conf_contour_hist_norm = BH3_hist_norm;
else
    conflict.sc_conf_contour = zeros(1,60);
    conflict.sc_conf_contour_hist = zeros(1,60);
    conflict.sc_conf_contour_hist_norm = zeros(1,60);
end


