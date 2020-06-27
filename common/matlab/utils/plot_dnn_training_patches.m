format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% select the data file

file_filter = {'*.bin','binary files (*.bin)'; '*.mat','MAT-files (*.mat)'; '*.txt','Text Files'; '*.*','All Files' };

startpath = 'D:\IUPUI\PhD\Results\';
[data_file, data_path] = uigetfile(file_filter, 'Select Cropper Stats File', startpath);
if(data_path == 0)
    return;
end

commandwindow;

%% get the data

fprintf('Reading in cropper statistics...\n');

[~,data_file_name,file_ext] = fileparts(data_file);

if(strcmp(file_ext,'.txt')==1)
    data_params = parse_input_parameters(fullfile(data_path, data_file));

    % convert the cells to a matrix
    data = zeros(numel(data_params),5);
    for idx=1:numel(data_params)
        data(idx,:) = str2double(data_params{idx});
    end

elseif(strcmp(file_ext,'.mat')==1)
    
    load(fullfile(data_path, data_file));
    
elseif(strcmp(file_ext,'.bin')==1)
    
    data = read_binary_cropper_stats(fullfile(data_path, data_file));
    %save(fullfile(data_path, strcat(data_file_name,'.mat')),'data');
else
    return;
end

save_path = data_path;  % 'D:\IUPUI\PhD\Images\dfd_dnn\';    

%% one more look at the distribution of actual depth map values selected

% select the data input file that contains the locations of the depth maps
file_filter = {'*.txt','Text Files'; '*.*','All Files' };

startpath = 'D:\IUPUI\DfD\';
[train_data_file, train_data_path] = uigetfile(file_filter, 'Select the training input file', startpath);
if(train_data_path == 0)
    return;
end

[params] = parse_input_parameters(fullfile(train_data_path, train_data_file));

data_directory = strcat('D:\IUPUI\',params{1}{1}); 
params(1) = [];

commandwindow;

%% run through and get some info about the training patches

img_count = numel(params);
hist_bins = 0:1:img_count;

h = histcounts(data(:,1),hist_bins);

min_hist = min(h);
mean_hist = mean(h);
max_hist = max(h);

fprintf('Minimum Image Selection: %4.2f\n', min_hist);
fprintf('Average Image Selection: %4.2f\n', mean_hist);
fprintf('Maximum Image Selection: %4.2f\n\n', max_hist);

%% Plot the histogram of the patch selection

figure(plot_num);
set(gcf,'position',([50,50,1200,500]),'color','w')
hold on
box on
grid on

b = bar(hist_bins(1:end-1),h,'b');

b(1).FaceColor = 'b';

% X-Axis
set(gca,'TickLabelInterpreter','none', 'fontweight', 'bold')
xlim([hist_bins(1)-1 hist_bins(end)]);
xticks([0:20:img_count-1]);
xlabel(strcat('Training Image Number'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
% plt_max = 0.9912;
% plt_min = 0.989;
% ylim([plt_min plt_max]);
% yticks([plt_min:0.0002:plt_max]);
% ytickformat('%1.4f')
ylabel('# of Training Patches From Image', 'fontweight', 'bold', 'FontSize', 13);

title('DfD-Net Training Patch Distribution','fontweight', 'bold', 'FontSize', 13);

ax = gca;
ax.Position = [0.07 0.11 0.9 0.82];

print(plot_num, '-dpng', fullfile(save_path, strcat('dfd_dnn_train_patch_dist.png')));

plot_num = plot_num + 1;

%% plot an example image

if(false)
    
    crop_w = 32;
    crop_h = 32;

    crops = data(data(:,1)==400,:);

    a = 1/(size(crops,1));

    img_h = crops(1,2);
    img_w = crops(1,3);

    figure(plot_num);
    set(gcf,'position',([50,50,800,600]),'color','w')
    hold on
    box on
    grid on

    for idx=1:size(crops,1)

        v = [crops(idx,4),crops(idx,5);
             crops(idx,4)+crop_w,crops(idx,5);
             crops(idx,4)+crop_w,crops(idx,5)+crop_h;
             crops(idx,4),crops(idx,5)+crop_h];

        patch('Faces',[1 2 3 4],'Vertices',v,'FaceColor','green','FaceAlpha',0.005, 'EdgeColor','green','LineWidth',0.5, 'LineStyle','-'); 

    end

    v = [0,0; img_w-1,0; img_w-1,img_h-1;0,img_h-1];
    patch('Faces',[1 2 3 4],'Vertices',v,'FaceColor','none','FaceAlpha',1, 'EdgeColor','k','LineWidth',1, 'LineStyle','-');

    xlim([0 img_w-1]);

    ylim([0 img_h-1]);

    axis off;

    ax = gca;
    ax.Position = [0.04 0.04 0.92 0.92];

    print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_dnn_train_patch_coverage.png')));

    plot_num = plot_num + 1;

end

%% create a 2-D surface of the points
jet_max = 110;
crop_w = 32;
crop_h = 32;

use_min = false;

if(use_min)
    [~, img_index] = min(h);
    plot_name = 'min';
else
    [~, img_index] = max(h); 
    plot_name = 'max';
end


crops = data(data(:,1)==(img_index-1),:);

% adjust for matlab indexing
crops(:,4) = crops(:,4) + 1;
crops(:,5) = crops(:,5) + 1;


img_h = crops(1,2);
img_w = crops(1,3);

patch = ones(crop_h, crop_w);

patch_map = zeros(img_h, img_w);
for idx=1:size(crops,1)

    x = crops(idx,4)+1;
    y = crops(idx,5)+1;
    
    % reverse y axis
    %patch_map(y:-1:y-crop_h+1,x:x+crop_w-1) = patch_map(y:-1:y-crop_h+1,x:x+crop_w-1) + patch;
    
    % standard y axis
    if((y+crop_h-1 < img_h) && (x+crop_w-1 < img_w))
        patch_map(y:1:y+crop_h-1,x:x+crop_w-1) = patch_map(y:1:y+crop_h-1,x:x+crop_w-1) + patch;
    end
end

figure(plot_num);
set(gcf,'position',([50,50,800,600]),'color','w')
hold on
box on
grid on

%image(patch_map);
imagesc(patch_map);
colormap(jet(jet_max));

% X-Axis
xlim([1 img_w]);
%xticks([0:20:img_w-1]);
xticklabels([]);

% Y-Axis
ylim([1 img_h]);
yticks([0:20:img_h-1]);
yticklabels([]);

%colormap(jet(jet_val));
cb = colorbar('fontweight','bold','FontSize', 13, 'Location', 'eastoutside');
%[cb.Ticks , cb.TickLabels] = calc_limits(0, nmae_plt_max, plt_step, '%1.2f');
cb.Label.String = '# of Pixels Chosen';
cb.Ticks = [0:10:jet_max];
cb.TickLabels = num2str(cb.Ticks');
cb.Limits = [0 jet_max];

ax = gca;
ax.Position = [0.01 0.02 0.87 0.96];

print(plot_num, '-dpng', fullfile(data_path, strcat('tp_hm_plt_', plot_name, num2str(img_index,'_%03d'), '.png')));

plot_num = plot_num + 1;


%% convert the patch heat map to an image

RGB = ind2rgb(patch_map,colormap(jet(jet_max))); 

save_path = 'D:\IUPUI\PhD\Images\test_v3';
mkdir(save_path);

%imwrite(RGB,fullfile(save_path, strcat('training_patch_heatmap_', num2str(img_index,'%03d'), '.png')));
imwrite(RGB,fullfile(data_path, strcat('tp_hm_', plot_name, num2str(img_index,'_%03d'), '.png')));


%% read in the depth map images

dm_img = cell(img_count,1);
gt_hist = zeros(256,1);

for idx=1:img_count
    
    dm_img{idx,1} = imread(fullfile(data_directory,params{idx}{3}));
    
    img_size = mod(size(dm_img{idx,1}),16);
    dm_img{idx,1} = dm_img{idx,1}(1:end-img_size(1),1:end-img_size(2));
    
    tmp = dm_img{idx,1}(:);
%     tmp = tmp(:);

    for jdx=1:256
        gt_hist(jdx,1) = gt_hist(jdx,1) + sum(tmp==(jdx-1));        
    end
%     for jdx=1:numel(tmp)
%         gt_hist(tmp(jdx)+1,1) = gt_hist(tmp(jdx)+1,1) + 1;
%     end
    
end
fprintf('Done reading in images!\n');

%% cycle through the data
tic;
dm_hist = zeros(256,1);
for idx=1:size(data,1)
    
    img_num = data(idx,1)+1;
    
    x = data(idx,4)+1;
    y = data(idx,5)+1;
      
    if((y+crop_h-1 < size(dm_img{img_num,1},1)) && (x+crop_w-1 < size(dm_img{img_num,1},2)))
        tmp = dm_img{img_num,1}(y:1:y+crop_h-1,x:x+crop_w-1);
        tmp = tmp(:);

%         for jdx=1:256
%             dm_hist(jdx,1) = dm_hist(jdx,1) + sum(tmp==(jdx-1));        
%         end
        for jdx=1:numel(tmp)
            dm_hist(tmp(jdx)+1,1) = dm_hist(tmp(jdx)+1,1) + 1;
        end
    end      
    
end
toc;
fprintf('Done computing depth map usage!\n');


%% plot the hist of the dm values selected
hist_bins = 0:1:256;

% calculate the ratio between the two
dm_ratio = dm_hist./gt_hist;
nn = ~isnan(dm_ratio);

fprintf('Min Ratio:     %2.2f\n', max(dm_ratio(nn)));
fprintf('Average Ratio: %2.2f\n', mean(dm_ratio(nn)));
fprintf('Max Ratio:     %2.2f\n', min(dm_ratio(nn)));



h1 = cat(2,dm_hist,zeros(256,1));
h2 = cat(2,zeros(256,1),gt_hist);
h3 = cat(2, gt_hist, dm_hist);

figure(plot_num)
set(gcf,'position',([100,100,1300,600]),'color','w')
hold on
box on
grid on

%yyaxis left
b1 = bar(hist_bins(1:end-1), h3);
% b1 = bar(hist_bins(1:end-1), dm_ratio);

set(gca,'fontweight','bold','FontSize',13,'yscale','log');


% yyaxis right
% b2 = bar(hist_bins(1:end-1),h2);
% set(gca,'fontweight','bold','FontSize',13);

% X-Axis
xlim([-2, 257]);
xticks([0:5:255]);
xtickangle(90);
xlabel(strcat('Depth Map Value Count'),'fontweight','bold')

% Y-Axis
ylabel('Depth Map Value Count','fontweight','bold');

b1(1).FaceColor = 'b';
b1(2).FaceColor = 'r';

title('Ground Truth/Training Distribution', 'fontweight','bold','FontSize',16);
%lgd = legend([b1(1), b1(2)],'Training Data','Testing Data', 'location','southoutside', 'orientation', 'horizontal');
lgd = legend([b1(1), b1(2)],'Ground Truth Data', 'Training Sample Data', 'location','southoutside', 'orientation', 'horizontal');

ax = gca;
ax.Position = [0.055 0.17 0.93 0.77];

print(plot_num, '-dpng', fullfile(data_path,'dm_combined.png'));

plot_num = plot_num + 1;

%% save the data in the form of the dist viewer mat file

save_path = 'D:\IUPUI\PhD\Images\datasets';

file_filter = {'*.mat','MAT-files (*.mat)'; '*.bin','binary files (*.bin)'; '*.*','All Files' };
[mat_file, mat_path] = uiputfile(file_filter, 'Select Mat File Save Location', data_path);

if(mat_path == 0)
    return;
end

train_hist = dm_hist'; %h1(:,1)';
test_hist = [];
name = 'Training Samples';
train_hist_sum = sum(dm_hist(:));       %size(data,1)*(32*32);
test_hist_sum = 0;

save(fullfile(mat_path, mat_file), 'train_hist', 'train_hist_sum', 'test_hist', 'test_hist_sum', 'name');
            

return;

%% new plot that overlays the actual distributions

h1 = cat(2,dm_hist,zeros(256,1));
h2 = cat(2,zeros(256,1),gt_hist);

figure(plot_num)
set(gcf,'position',([100,100,1300,600]),'color','w')
hold on
box on
grid on

yyaxis left
b1 = bar(hist_bins(1:end-1), h1);
set(gca,'fontweight','bold','FontSize',13);


yyaxis right
b2 = bar(hist_bins(1:end-1),h2);
set(gca,'fontweight','bold','FontSize',13);

% X-Axis
xlim([-2, 257]);
xticks([0:5:255]);
xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ylabel('Depth Map Ratio','fontweight','bold');

b1(1).FaceColor = 'b';
b2(2).FaceColor = 'g';

title('Ground Truth/Training Distribution Ratio', 'fontweight','bold','FontSize',16);
lgd = legend([b1(1), b2(2)],'Training Data','Testing Data', 'location','southoutside', 'orientation', 'horizontal');

ax = gca;
ax.Position = [0.05 0.15 0.90 0.80];

%print(plot_num, '-dpng', fullfile(save_dir,'dm_combined.png'));

plot_num = plot_num + 1;
