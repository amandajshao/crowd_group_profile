function fun_grVisual(visual_data, path_img, frames, t0, t1, color_ind, t_all)

% FUN_GRVISUAL: Group visualization
%               Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% visualization
if nargin == 4
    % there is no coherent trk at time t0, only show the image
    f1=figure(1);
    curFrame = imread([path_img frames(t0).name]);
    imshow(curFrame)
    title(['Frame No.' num2str(t0)])
elseif nargin == 6
    f1=figure(1);
    tt = 0;
    for t = t0 : t1
        curAllX = [];
        curAllIndex = [];
        tt = tt + 1;
        curFrame = imread([path_img frames(t).name]);
        imshow(curFrame)
        title(['Frame No.' num2str(t)])
        hold on
        % search trk of each time
        num_temp = 0;
        nCluster = sum(~cellfun(@isempty, visual_data(:,tt)));
        for i = 1 : nCluster
            for j = 1 : length(visual_data{i,tt})
                if ~isempty(find(visual_data{i,tt}(j,1).t==t, 1))
                    num_temp = num_temp + 1;
                    % store each trk's coordinate at this time
                    curAllX(1,num_temp) = visual_data{i,tt}(j,1).x(visual_data{i,tt}(j,1).t==t);
                    curAllX(2,num_temp) = visual_data{i,tt}(j,1).y(visual_data{i,tt}(j,1).t==t);
                    % store the group index of the above each trk
                    curAllIndex(1,num_temp) = i;
                end
            end
            if ~isempty(curAllIndex)
                cr = fun_colorrand(color_ind{:,tt}(i,1));
                plot(curAllX(1,curAllIndex==i), curAllX(2,curAllIndex==i), ...
                    'o', 'MarkerSize', 5, 'MarkerEdgeColor', cr, ...
                    'MarkerFaceColor', cr)
                hold on
            end
        end
        drawnow
        hold off
    end
elseif nargin == 7  
    f1=figure(1);
    for t = 1 : length(t_all)
        if t_all{1,t} ~= 0
            curAllX = [];
            curAllIndex = [];
            curFrame = imread([path_img frames(t).name]);
            imshow(curFrame)
            title(['Frame No.' num2str(t)])
            hold on
            % search trk of each time
            num_temp = 0;
            nCluster = sum(~cellfun(@isempty, visual_data(:,t)));
            for i = 1 : nCluster
                for j = 1 : length(visual_data{i,t})
                    if ~isempty(find(visual_data{i,t}(j,1).t==t, 1))
                        num_temp = num_temp + 1;
                        % store each trk's coordinate at this time
                        curAllX(1,num_temp) = visual_data{i,t}(j,1).x(visual_data{i,t}(j,1).t==t);
                        curAllX(2,num_temp) = visual_data{i,t}(j,1).y(visual_data{i,t}(j,1).t==t);
                        % store the group index of the above each trk
                        curAllIndex(1,num_temp) = i;
                    end
                end
                if ~isempty(curAllIndex)
                    cr = fun_colorrand(color_ind{:,t}(i,1));
                    plot(curAllX(1,curAllIndex==i), curAllX(2,curAllIndex==i), ...
                        'o', 'MarkerSize', 5, 'MarkerEdgeColor', cr, ...
                        'MarkerFaceColor', cr)
                    hold on
                end
            end
            drawnow
            hold off
        end
    end
end


end