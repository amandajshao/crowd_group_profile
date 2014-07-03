function [points_trans, points_trans_normalize, density_pattern] = fun_calibration(groupPoint, ...
    smooth_size, smooth_sigma, pattern_scale, gv_flag)

% FUN_CALIBRATION: Find affine matrix to normalize all groups into one axis system
%                  Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% calibration
meanV = groupPoint(:,4:5);
invalid_idx = sum(abs(meanV),2)==0;
meanV(invalid_idx,:) = [];
meanV = mean(meanV);

% scale vector
s = norm(meanV);

% rotation vector
cos_theta = meanV(2)/s; sin_theta = meanV(1)/s;
R = [cos_theta, -sin_theta; sin_theta, cos_theta];
% theta = acos(dot(meanV, vNormal)./(norm(meanV)));
% theta = atan(meanV(1)/meanV(2));
% R = ([cos(theta), -sin(theta); sin(theta), cos(theta)]);

% translation vector
t = mean(groupPoint(:,2:3));  

% transform axis -- scale + rotation + translation
points = groupPoint(:,2:3);
points_trans = (bsxfun(@minus, points, t))';
if gv_flag == 1 % 1 for group; 0 for video
    points_trans = R*points_trans/s;
end

%% normalize
minXY = min(points_trans, [], 2); maxXY = max(points_trans,[],2);
lengthXY = maxXY - minXY;
lengthXY = max(lengthXY);
points_trans_normalize = bsxfun(@minus, points_trans, (minXY + maxXY)/2);
points_trans_normalize = bsxfun(@rdivide, points_trans_normalize, lengthXY);
% points_trans_normalize = bsxfun(@minus, points_trans_normalize, 0.5);

%% gaussian smooth
points_trans_smooth = round(points_trans_normalize*pattern_scale) + smooth_size + round(pattern_scale/2);
density_pattern = zeros(pattern_scale+smooth_size*2+1);
density_pattern(sub2ind([pattern_scale+smooth_size*2+1, pattern_scale+smooth_size*2+1], points_trans_smooth(2,:)+1, points_trans_smooth(1,:)+1)) = 1;

H = fspecial('gaussian', smooth_size, smooth_sigma);
density_pattern = imfilter(density_pattern, H);

