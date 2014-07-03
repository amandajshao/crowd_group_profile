function draw_cm(mat,tick,num_class)
%%
%  Matlab code for visualization of confusion matrix;
%  Parameters£ºmat: confusion matrix;
%              tick: name of each class, e.g. 'class_1' 'class_2'...
%              num_class: number of class
%
%  Author£º Page( Ø§×Ó)  
%           Blog: www.shamoxia.com;  
%           QQ:379115886;  
%           Email: peegeelee@gmail.com
%%
set(gcf, 'color','w');
imagesc(1:num_class,1:num_class,mat);            %# in color
colormap(flipud(gray));  %# for gray; black for large value.
% colormap(jet);
colorbar; caxis([0 1]);

textStrings = num2str(mat(:),'%0.2f');  
textStrings = strtrim(cellstr(textStrings)); 
[x,y] = meshgrid(1:num_class); 
hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center', ...
    'FontSize', 20, 'FontWeight','bold', 'FontName', 'Times New Roman');%
midValue = mean(get(gca,'CLim')); 
% midValue = 0;
textColors = repmat(mat(:) >= midValue,1,3); 
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(gca,'xticklabel',tick,'XAxisLocation','top', ...
    'FontSize', 20, 'FontName', 'Times New Roman');%'FontWeight','bold', 
set(gca, 'XTick', 1:num_class, 'YTick', 1:num_class, ...
    'FontSize', 20, 'FontName', 'Times New Roman');%'FontWeight','bold', 
set(gca,'yticklabel',tick, 'FontSize', 20, 'FontName', 'Times New Roman');%'FontWeight','bold', 
rotateXLabels(gca, 45 );% rotate the x tick 315



