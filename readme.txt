Source code for Scene Independent Group Profiling in Crowd.

%======================= Main functions ========================%
1, main_gr.m
If you want to run code for one video, from group detection to group descriptor, please run this main function.
2, main_video.m
If you want to run code for multi-videos, from group detection to group descriptors, and then video classification, please run this main function.
3, main_video_class.m
This is for directly checking the classification results with saved middle results.

%====================== Functions for each part ================%
1, main_gr_detection.m 
Main function of group detection for one specific group and check the result
2, main_video_detection.m
Main function of group detection for multi-videos
3, main_gr_collectiveness.m/ main_gr_stability.m/ main_gr_uniform.m/ main_gr_conflict.m/ main_size.m
Core function of group descriptor for one video.
4, main_video_collectiveness.m/ main_video_stability.m/ main_video_uniform.m/ main_video_conflict.m+ main_video_BoW.m/ main_video_size.m
Core function of group descriptors for multi-videos ("Collectiveness","Stability","Uniform","Conflict", "GroupSize")


%======================= Relevant functions and instructions ====================%
1, cohGr_det.m
Core function of group detection algorithm. It follows algorithm of 4 steps in the paper:
step1: generate coherent filtering clusters [fun_CF.m]
step2: identify anchor tracklet [fun_anchorTrk.m]
step3: discover seeding tracklets [fun_groupInit.m] for learning CT prior [fun_kf.m]
step4: group refinement [fun_groupRef.m]

2, Details of group detection can be found in "CoheGroupDetect" folder

3, Details of group descriptors can be found in "descriptor" folder

4, Some utility functions are included in "util" folder.

5, No mex need before running code, just add "CoheGroupDetect", "descriptor", and "util" folder as well as their subfolders to the path.

6, One example image sequence is provided in "1_8_groupSplit-festivalwalk_1_2-1". You can first run main_gr_detection to get group detection results, the collective transition prior A and other elements for group descirptor computation. 

7, In this version, I put 5 descriptors separatively in 5 m-files. It is easy to join them together by yourself.

8, There are two ground truth file and one video information file used in the code, they are
- gt_labels.mat: groundtruth for video labels
- gt_videoName_bool.mat: for leave-one-out testing
- video_info_t0.xls: video information [video name,the video length, video size, video source, start time for group descriptor computation, end time for group descriptor computation]


-----------------------------------------------------------------------
July. 3, 2014, Jing Shao

If you use our code, please cite our papers:
Jing Shao, Chen Change Loy, and Xiaogang Wang. Scene-Independent Group Profiling in Crowd. CVPR, 2014.
Jing Shao, Chen Change Loy, and Xiaogang Wang. Learning Scene-Independent Group Descriptors for Crowd Understanding. TCSVT, 2016.  

