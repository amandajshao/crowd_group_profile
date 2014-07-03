function [CF_init, len_CF] = fun_CFUpd(CF_init, t, id_update)

% FUN_CFUPD: Update CF_init once a group has been detected
%            Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

if nargin == 3
    if ~isempty(CF_init{1,t})
        for i = 1 : length(CF_init{1,t}.gr)
            if ~isempty(intersect(CF_init{1,t}.gr{i,1}(:,1), id_update))
                [~,id_cf_temp,~] = intersect(CF_init{1,t}.gr{i,1}(:,1), id_update);
                CF_init{1,t}.gr{i,1}(id_cf_temp,1) = 0;
                CF_init{1,t}.gr{i,1}(id_cf_temp,2) = -1;
            end
        end
    end
    len_CF = zeros(length(CF_init{1,t}.gr),1);
    for i = 1 : length(CF_init{1,t}.gr)
        len_CF(i,1) = len_CF(i,1) + length(find(CF_init{1,t}.gr{i,1}(:,1)~=0));
    end
elseif nargin == 2
    len_CF = zeros(length(CF_init{1,t}.gr),1);
    for i = 1 : length(CF_init{1,t}.gr)
        len_CF(i,1) = len_CF(i,1) + length(find(CF_init{1,t}.gr{i,1}(:,1)~=0));
    end
end

end

