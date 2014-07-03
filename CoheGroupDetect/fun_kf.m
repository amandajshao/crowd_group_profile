function [A, x_train] = fun_kf(group_seed, max_iter, t, t_interval)

% FUN_KF: Kalman filtering
%         Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% training data
x_train = [];
for i = 1 : length(group_seed)
    x_train{i,1}(1,:) = group_seed(i,1).x(group_seed(i,1).t>=t & ...
        group_seed(i,1).t<=t+t_interval-1);
    x_train{i,1}(2,:) = group_seed(i,1).y(group_seed(i,1).t>=t & ...
        group_seed(i,1).t<=t+t_interval-1);
    x_train{i,1}(3,:) = ones(1,length(x_train{i,1}(1,:)));
end
seq_num = length(x_train);

%% parameters and initialization
R = diag([1,1,0]);
Q = R;
R = 0.1*R;
H = eye(3);

% use least square to evaluate the initial F
B1 = zeros(3, seq_num);
B2 = B1;
for i = 1 : seq_num
    B1(:,i) = [x_train{i,1}(1,1); x_train{i,1}(2,1); 1];
    B2(:,i) = [x_train{i,1}(1,2); x_train{i,1}(2,2); 1];
end
F = (pinv(B1*B1') * (B1*B2'))';

initx = sum(B1,2)/seq_num;
initx = initx(:,1);

initV = 1*eye(3);

%% Kalman filter learning
[A, ~, ~, ~, ~, ~, ~] = learn_kalman(x_train, F, H, Q, R, initx, initV, max_iter);

end