% training process
% X(t+1) = F X(t) + noise(Q)
% Y(t) =  H Y(t) + noise(R)


function [xss, x_train, F2] = ...
    em_training(trkTrack, CF_trk, jpgFolder, imTrkIncom, trange, max_iter, trk_num)

% preprocess CF_trk
[CF_trk_new, tstart, tend, x_train, x_obs_part, P_part] = pre_CFtrk(CF_trk, trkTrack, trange);
seq_num = length(CF_trk_new);
ss = seq_num*3; % state size
os = seq_num*3; % observation size
xss = [];
x_train_n = zeros(3,size(x_train,2),seq_num);
% F_all = zeros(3,3,seq_num);
% Q_all = zeros(3,3,seq_num);
% initx_all = zeros(3,seq_num);
% initV_all = zeros(3,3,seq_num);

if length(CF_trk_new) < 2
%     savefile = ['D:\jshao\dataset_mine\ICCV_video_experiment\4_1_NoPrediction\', ...
%         num2str(trk_num), '.jpg'];
    I = imread([jpgFolder, sprintf('%06d',trkTrack.t(end)), '.jpg']);
    f5 = figure(5); 
%     set(f5, 'visible', 'off');
    imshow(I); hold on;
    plot(trkTrack.x, trkTrack.y, 'o', 'MarkerSize', 1.5, ...
        'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')
%     saveas(gcf, savefile)

    % take a mark of whether there is a group for the trkTrack or not
    x_train = [];
    F2 = [];
else
    %% test show coherent tracks used for predicting
    f5 = figure(5); 
%     set(f5, 'visible', 'off');
    imshow(imTrkIncom); hold on;
    for i = 1 : seq_num
        plot(x_train(3*i-2,:), x_train(3*i-1,:), '.g')
        hold on
    end
    hold off
    
    %% define parameters
    for n = 1 : seq_num
        x_train_n(:,:,n) = x_train(3*n-2:3*n, :);
    end
    R = diag([1,1,0]);
    Q = R;
    R = 0.1*R;
    H = eye(3);
    
    % use least square to evaluate the initial F
    B1 = zeros(3, seq_num);
    B2 = B1;
    B1(1,:) = x_train(1:3:end,1)';
    B1(2,:) = x_train(2:3:end,1)';
    B1(3,:) = ones(1,seq_num);
    B2(1,:) = x_train(1:3:end,2)';
    B2(2,:) = x_train(2:3:end,2)';
    B2(3,:) = ones(1,seq_num);
    F = ((B1*B1')\(B1*B2'))';
    
    initx = sum(x_train_n,3)/seq_num;
    initx = initx(:,1);
    
    initV = 1*eye(3);
    
    %% Kalman filter learning
    [F2, H2, Q2, R2, initx2, initV2, LL] = ...
        learn_kalman(x_train_n, F, H, Q, R, initx, initV, max_iter);
    % force Q into a certain form with 0
    % Here, Q2(3*i,:)=0 && Q2(:,3*i)=0 or Q2(3*i,3*i)=0 are same
%     for i = 1 : ss/3
%         Q2(3*i,3*i) = 0;
%         % Q2(3*i,:) = 0;
%         % Q2(:,3*i) = 0;
%     end
    
    %% Kalman smoother on the incomplete trk
    [xsmooth, Vsmooth] = kalman_smoother(x_train_n(:,:,1), F2, H2, Q2, R2, initx2, initV2);
    % savefold = [jpgFolder, num2str(trk_num)];
    % mkdir(savefold);
    xss = [];
    xnew = xsmooth(:,end);
        
    %% show the result    
    frames = dir([jpgFolder '*.jpg']);
    nFrame = length(frames);
    
    f6 = figure(6); 
%     set(f6, 'visible', 'off');
    if trkTrack.t(1) ~= tstart
        for t_temp = trkTrack.t(1) : tstart-1
            I = imread([jpgFolder, sprintf('%06d',t_temp), '.jpg']);
            imshow(I); hold on;
            
            plot(trkTrack.x(1:find(trkTrack.t==t_temp)), ...
                trkTrack.y(1:find(trkTrack.t==t_temp)), ...
                'o', 'MarkerSize', 1.5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
            hold on; drawnow; hold off;
            % savefile = [savefold, '\', num2str(t_temp-trkTrack.t(1)+1), '.jpg'];
            % saveas(gcf, savefile)
        end
    end
    
    % options = optimset('Algorithm', 'Levenberg-Marquardt');
    
    for t_temp = tstart : min(tend+trange, nFrame-1)
        I = imread([jpgFolder, sprintf('%06d',t_temp), '.jpg']);
        imshow(I); hold on;
        % plot(trkTrack.x, trkTrack.y, 'o', 'MarkerSize', 1.5, ...
        %     'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g')
        % hold on
        if trkTrack.t(1) ~= tstart
            plot(trkTrack.x(1:find(trkTrack.t==(tstart-1))), ...
                trkTrack.y(1:find(trkTrack.t==(tstart-1))), ...
                'o', 'MarkerSize', 1.5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        end
        
        % show the smooth result
        if t_temp <= tend
            xs = xsmooth(1:3, t_temp-tstart+1);
            xss = [xss, xs];
            plot(xss(1,:), xss(2,:), 'o', 'MarkerSize', 1.5, ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
            hold on; drawnow; hold off;
        % show the predicted result
        else
            % no long observation help to predict
%             if isempty(x_obs_part) 
                xnew = F2*xnew; % x_t+1 = F*x_t
            % x_t+1 predicted by combining partial observations
%             else
%                 xnew = pinv(pinv(Q2)+P_part'*pinv(P_part*R2*P_part')*P_part) * ...
%                     (P_part'*pinv(P_part*R2*P_part')*x_obs_part(:,t_temp-tend) + ...
%                     pinv(Q2)*F2*xnew);
% %                 xnew = (inv(Q)+P_part'/(P_part*R*P_part')*P_part) \ ...
% %                     (P_part'/(P_part*R*P_part')*x_obs_part(:,t_temp-tend) + Q\F2*xnew)
% %                 xnew = lsqnonlin(@(x) myfun(x, P_part, R2, F2, Q2, x_obs_part(:,t_temp-tend), xnew), xnew, [], [], options);
%             end
            xs = xnew(1:3, :);
            xss = [xss, xs];
            plot(xss(1,1:tend-tstart), xss(2,1:tend-tstart), 'o', ...
                'MarkerSize', 1.5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
            hold on
            plot(xss(1,tend-tstart:t_temp-tstart), xss(2,tend-tstart:t_temp-tstart), ...
                'o', 'MarkerSize', 1.5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')
            drawnow; hold off;
        end
        
        % savefile = [savefold, '\', num2str(t_temp-trkTrack.t(1)+1), '.jpg'];
        % saveas(gcf, savefile)
    end
    
    % store the complete track
    f7 = figure(7); 
%     set(f7, 'visible', 'off');
    sp1 = subplot(1,2,1);
    pos1 = get(sp1, 'pos');
    pos1(1) = pos1(1) - 0.1;
    pos1(2) = pos1(2) - 0.1;
    pos1(3) = pos1(3) + 0.13;
    pos1(4) = pos1(4) + 0.13;
    set(sp1, 'pos', pos1)
    I = imread([jpgFolder, sprintf('%06d',trkTrack.t(end)), '.jpg']);
    imshow(I); hold on;
    plot(trkTrack.x,  trkTrack.y, 'o', 'MarkerSize', 1.5, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
    
    if trkTrack.t(1) ~= tstart
        temp = [trkTrack.x(1:find(trkTrack.t==(tstart-1)))'; ...
            trkTrack.y(1:find(trkTrack.t==(tstart-1)))'; ...
            ones(1,find(trkTrack.t==(tstart-1)))];
        xss = [temp, xss];
    end
    
    sp2 = subplot(1,2,2);
    pos2 = get(sp2, 'pos');
    pos2(1) = pos2(1) - 0.1;
    pos2(2) = pos2(2) - 0.1;
    pos2(3) = pos2(3) + 0.13;
    pos2(4) = pos2(4) + 0.13;
    set(sp2, 'pos', pos2)
    I = imread([jpgFolder, sprintf('%06d',t_temp), '.jpg']);
    imshow(I); hold on;
    plot(xss(1,:), xss(2,:), 'o', 'MarkerSize', 1.5, ...
        'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y')

%     savefile = ['C:\Users\lsheng\Dropbox\sj\4_1_ComTrkPred_Separate\', ...
%         num2str(trk_num), '.jpg'];
%     saveas(gcf, savefile)

    hold off
end



