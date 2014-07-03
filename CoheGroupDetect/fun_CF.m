function [CF_init,trks] = fun_CF(trks, path_img, t_begin, t_end, paraSet)

% FUN_CF: Main function for coherent filtering
%         Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% parameter setting
d = paraSet.d;
K = paraSet.K;
lamda = paraSet.lamda;
[trkTime, ~, ~] = fun_trkInfo(trks);

CF_init = [];
f1=figure(1);
set(f1,'visible','off');
for i = t_begin : t_end    
    fprintf('CF frame %d.\n', i);
    % coherent filtering clustering
    curIndex = (find(trkTime(1,:)<=(i) & trkTime(2,:)>=(i)))';    
    includedSet = trks(1,curIndex);
    [curAllX,clusterIndex] = CoherentFilter(includedSet,i,d,K,lamda); 
    nCluster = max(clusterIndex);

    % store clusters and show result
    if ~isempty(curAllX)
        for j = 1 : length(curAllX(3,:))
            curAllX(3,j) = curIndex(curAllX(3,j));
        end
        for j = 1 : nCluster
            % store index of each group in every frame, denoted by "CF_init"
            CF_init{1,i-t_begin+1}.gr{j,1}(:,1) = curAllX(3,clusterIndex==j)';
            CF_init{1,i-t_begin+1}.gr{j,1}(:,2) = zeros(length(find(clusterIndex==j)),1);% record the trk can be anchor or not
            for jj = 1 : length(CF_init{1,i-t_begin+1}.gr{j,1})
                trks(CF_init{1,i-t_begin+1}.gr{j,1}(jj,1)).num = ...
                    trks(CF_init{1,i-t_begin+1}.gr{j,1}(jj,1)).num + 1;
            end
        end
    end
end

end






