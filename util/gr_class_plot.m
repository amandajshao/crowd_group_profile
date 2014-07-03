function [accuracy_mean, C, confusion_matrix] = gr_class_plot(label, result_fig_path, path)
% label: [test_id, gt_label, test_label]

accuracy_mean = []; C = []; confusion_matrix = [];

%% do confusion matrix and plot -- version 1
[C,order] = confusionmat(label(:,2),label(:,3));
confusion_matrix = C./repmat(sum(C,2),1,size(C,1));

accuracy_mean = mean(diag(confusion_matrix));

% figure(2); 
% set(gcf, 'color','w');
% imagesc(confusion_matrix); colorbar; caxis([0,1])
% 
% saveFile = [result_fig_path, 'confusionMatrix_v1.pdf'];
% set(gcf,'paperpositionmode','auto');
% saveas(gcf, saveFile, 'pdf');
% 
% % saveFile = [result_fig_path, 'confusionMatrix.eps'];
% % saveas(gcf, saveFile, 'psc2')

%% do confusion matrix and plot -- version 2
% class_num = unique(label(:,2));
% class_result = [];
% 
% for class_n = 1 : max(class_num)
%     class_gt_id = find(label(:,2) == class_n);
%     class_gt = label(class_gt_id, 2);
%     num_in_class(class_n,1) = length(class_gt);
%     class_result = [class_result; label(class_gt_id, 3)];
% %     name_class{class_n,1} = ['class',num2str(class_n)];
% end
% name_class = [{'Gas'}, {'Solid'}, {'Pure Liquid'}, {'Impure Liquid'}];
% 
% confusion_matrix = compute_confusion_matrix(class_result,num_in_class,name_class)
% accuracy_mean = mean(diag(confusion_matrix));
% 
% % saveFile = [result_fig_path, 'confusionMatrix_v2.pdf'];
% % set(gcf,'paperpositionmode','auto');
% % saveas(gcf, saveFile, 'pdf');
% % 
% % % saveFile = [result_fig_path, 'confusionMatrix_v2.eps'];
% % % set(gcf,'paperpositionmode','auto');
% % % print('-depsc', saveFile);
% % % saveas(gcf, saveFile, 'psc2')


%% show combination classification result (group jpg)
% % gr_path = [path, 'first_frame_gr(withGrDet)\'];
% % gr_path = [path, 'first_frame_gr(withGrDet_arrow)\'];
% % gr_file = dir([gr_path, '*.eps']);
% gr_path = [path, 'first_frame_gr(withGrFric)\'];
% gr_file = dir([gr_path, '*.jpg']);
% 
% class_num = unique(label(:,2));
% 
% for class_n = 1 : max(class_num)
%     class_id = find(label(:,2) == class_n);
% %     class_id = find(label(:,3) == class_n);
%     
%     class_gr_id = label(class_id, 1);
%     class_gt = label(class_id, 2);
%     class_result = label(class_id, 3);
%     
%     save_path = [result_fig_path, 'class_', num2str(class_n), '\'];
%     mkdir(save_path)
%     
%     for n = 1 : length(class_gr_id)
%         copyfile([gr_path,gr_file(class_gr_id(n)).name], ...
%             [save_path,'fric_',num2str(class_result(n)),'_',gr_file(class_gr_id(n)).name])
%     end
%     
% end








