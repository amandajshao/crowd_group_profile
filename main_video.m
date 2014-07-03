%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% main function of group detection %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% Please make sure add three folders and their subfolders to the path before running code
% They are: CoheGroupDetect, descriptor, util
% CoheGroupDetect contains group detection code, except main function
% descriptor contains descriptor code, except main function
% util contains some utilized code

%% multi-video group detection
main_video_detection; % it will save mid-result for group-descriptor in the folder '.\result_groupDet_new\'

%% multi-video group descriptors
main_video_collectiveness; 
main_video_stability;
main_video_uniform;
main_video_conflict;
main_video_BoW;
main_video_size;

%% video classification
main_video_class_v2; 



