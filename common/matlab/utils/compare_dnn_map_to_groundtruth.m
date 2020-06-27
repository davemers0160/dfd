format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% select the ground truth depth map and the resulting depth maps from the dnn

file_filter = {'*.png','Image Files';'*.*','All Files' };

data_path = 'D:\IUPUI\Test_Data\';
[gt_file, gt_path] = uigetfile(file_filter, 'Select Ground Truth Depth Map File', data_path, 'MultiSelect', 'off');
if(gt_path == 0)
    return;
end

data_path = 'D:\IUPUI\Test_Data\';
[img_file, img_path] = uigetfile(file_filter, 'Select In-Focus Image File', gt_path, 'MultiSelect', 'off');
if(img_path == 0)
    return;
end

data_path = 'D:\IUPUI\Test_Data\';
[dimg_file, dimg_path] = uigetfile(file_filter, 'Select Out-of-Focus Image File', img_path, 'MultiSelect', 'off');
if(dimg_path == 0)
    return;
end

data_path = 'D:\IUPUI\PhD\results';
[dnn_dm_file, dnn_dm_path] = uigetfile(file_filter, 'Select DNN Depth Map File', data_path, 'MultiSelect', 'off');
if(dnn_dm_path == 0)
    return;
end

save_data = false;

save_path = uigetdir(dnn_dm_path, 'Select Save Folder');
if(save_path == 0)
    save_data = false;
end

%% load in the data

% this is the maximum depthmap value possible
max_depthmap_value = 40;

% read in the dnn files first
dnn_img = double(imread(fullfile(dnn_dm_path,dnn_dm_file)));
dnn_sz = size(dnn_img);

% read in the in-focus image
img = double(imread(fullfile(img_path,img_file)));
img = img(1:dnn_sz(1), 1:dnn_sz(2),:);

% read in the out-of-focus image
dimg = double(imread(fullfile(dimg_path,dimg_file)));
dimg = dimg(1:dnn_sz(1), 1:dnn_sz(2),:);

% read in the ground truth
gt_img = double(imread(fullfile(gt_path, gt_file)));
gt_img = gt_img(1:dnn_sz(1), 1:dnn_sz(2));


%% try this instead

[gt_error_value, combined_hist, cm] = get_depthmap_errors(gt_img, dnn_img, max_depthmap_value);


%% get the dnn depth map distribution

dm_tmp = dnn_img(:);
gt_tmp = gt_img(:);

for idx=1:max_depthmap_value+1
    dm_hist(idx,1) = sum(dm_tmp==(idx-1));   
    gt_hist(idx,1) = sum(gt_tmp==(idx-1));
end
combined_hist = cat(2,gt_hist,dm_hist);

min_bin = 0;
max_bin = numel(gt_hist);

plot_step = 1;
hist_bin_step = 1;
hist_bins = min_bin:hist_bin_step:max_bin;

x_lim = [0,max_depthmap_value+1];
x = [min_bin:hist_bin_step:(max_bin-1)];

y_m = ceil(log10(max(combined_hist(:))));
y_max = 10^y_m;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
b = bar(x, combined_hist);
set(gca,'fontweight','bold','FontSize',13, 'yscale', 'log');

% X-Axis
xlim(x_lim);
%xticks([(x_lim(1)+2):plot_step:(x_lim(2)-2)]);
xticks([0:plot_step:max_bin-1]);
xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ylabel('Depth Map Value Count','fontweight','bold')
ylim([1 y_max]);
%ytickformat('%1.2f');

b(1).FaceColor = 'b';
b(2).FaceColor = 'r';

%title('Depth Map Distribution Comparison', 'fontweight','bold','FontSize',16);
lgd = legend('Ground Truth', 'DNN Depth Map', 'location','southoutside', 'orientation', 'horizontal');
ax = gca;
ax.Position = [0.07 0.17 0.91 0.80];

plot_num = plot_num + 1;


%% create a mask where the dnn depthmap is xx% from the ground truth

gt_diff = abs(gt_img - dnn_img);
a = gt_diff/(max_depthmap_value);

if(true)
figure(plot_num);
set(gcf,'position',([50,50,900,600]),'color','w', 'Name','Ground Truth Image');
image(gt_img);
colormap(jet(max_depthmap_value+1));
axis off;
plot_num = plot_num + 1;

figure(plot_num);
set(gcf,'position',([50,50,900,600]),'color','w', 'Name','DNN Depth Map');
image(dnn_img);
colormap(jet(max_depthmap_value+1));
axis off;
plot_num = plot_num + 1;

figure(plot_num);
set(gcf,'position',([50,50,900,600]),'color','w', 'Name','GT-DM Diff');
image(gt_diff);
colormap(jet(max_depthmap_value+1));
axis off;
plot_num = plot_num + 1;

figure(plot_num);
set(gcf,'position',([50,50,900,600]),'color','w', 'Name','Image Error Overlay');
image(uint8(img), 'alphadata', 1-a)
axis off;
plot_num = plot_num + 1;

% imwrite(uint8(img), fullfile(dnn_dm_path,'img_a.png'), 'png', 'Alpha', 1-a);
end

%%
%dm_rng = [55:1:58];   % Midd1 & Plastic overlap
%dm_rng = [21, 22, 23, 24, 25, 36, 37, 38]; % Midd2 overlap
dm_rng = [4:1:11];   % rw4 - k00, k04... overlap
%dm_rng = [154:1:159]; %rw4 - k47

gt_mask = zeros(dnn_sz);

for idx=1:numel(dm_rng)
    gt_mask = gt_mask | (gt_img == dm_rng(idx));
end
fprintf('Percent of mask values in depth map: %2.4f%%\n', (sum(gt_mask(:))*100)/numel(gt_mask(:)));


% figure(plot_num); 
% set(gcf,'position',([50,50,900,600]),'color','w', 'Name','Ground Truth Mask');
% image(gt_mask); 
% colormap(gray(2)); 
% axis off;
% plot_num = plot_num + 1;


%% combine the six images - no color

% make them the following:
% first row:  GT, DM, |GT-DM|
% second row: Mask, Mask/error, Image/error


    
% img_space_h = 20;     % number of pixels to put between the images
% img_size = [200,200];
% 
% pad = 255*ones(img_size(1), img_space_h);
% 
% gt_sm = imresize(gt_img, img_size, 'nearest');
% dnn_img_sm = imresize(dnn_img, img_size, 'nearest');
% gt_mask_sm = imresize(gt_mask*255, img_size, 'nearest');
% gt_diff_sm = imresize(gt_diff, img_size, 'nearest');
% 
% combined_img = cat(2, gt_sm, pad, dnn_img_sm, pad, gt_mask_sm, pad, gt_diff_sm);
% 
% if(false)
%     
% figure(plot_num); 
% set(gcf,'position',([50,50,1200,250]),'color','w', 'Name','Combined');
% image(combined_img); 
% colormap(gray(max_depthmap_value+1)); 
% axis off;
% plot_num = plot_num + 1;
% 
% %imwrite(uint8(combined_img),fullfile(dnn_dm_path,'combined_dm_error.png'));
% 
% end

%% make a color version of the six images

% make them the following:
% first row:  GT, DM, |GT-DM|
% second row: Mask, Mask/error, Image/error

img_space_h = 15;     % number of pixels to put between the images
img_size = [200,200];

pad = 255*ones(img_size(1), img_space_h, 3);

% grpund truth
gt_sm = imresize(gt_img, img_size, 'nearest');
%gt_sm = cat(3, gt_sm, gt_sm, gt_sm);
gt_sm = 255*ind2rgb(gt_sm, jet(max_depthmap_value+1));

% dnn depth map
dnn_img_sm = imresize(dnn_img, img_size, 'nearest');
%dnn_img_sm = cat(3, dnn_img_sm, dnn_img_sm, dnn_img_sm);
dnn_img_sm = 255*ind2rgb(dnn_img_sm, jet(max_depthmap_value+1));

% error
gt_diff_sm = imresize(gt_diff, img_size, 'nearest');
gt_diff_sm3 = cat(3, gt_diff_sm, gt_diff_sm, gt_diff_sm);
%gt_diff_sm = 255*ind2rgb(gt_diff_sm, jet(max_depthmap_value));
gt_diff_sm = (max_depthmap_value)*ind2rgb(gt_diff_sm, jet(max_depthmap_value+1));

% mask
gt_mask_sm = imresize(gt_mask*255, img_size, 'nearest');
gt_mask_sm = cat(3, gt_mask_sm, gt_mask_sm, gt_mask_sm);

% mask/error overlay
gt_mask_over = gt_diff_sm;
gt_mask_over(:,:,1) = max(gt_mask_over(:,:,1), gt_mask_sm(:,:,1));

% in-focus image/error overlay
img_sm = imresize(img, img_size, 'bicubic');
img_sm_over = (img_sm*0.5)+(gt_diff_sm3*0.5);

% out-of-focus image
dimg_sm = imresize(dimg, img_size, 'bicubic');

img_row1 = [img_sm pad dimg_sm pad gt_sm pad dnn_img_sm];
img_row2 = [gt_diff_sm pad gt_mask_sm pad gt_mask_over pad img_sm_over];

%% generate a plot that does not save the image, but saves the plot instead; needed for colorbar


figure(plot_num);
set(gcf,'position',([50,50,1300,300]),'color','w')

subplot(1,4,1);
imagesc(uint8(img_sm));
axis off;

ax = gca;
ax.Position = [0.00 0.00 0.24 1.00];

subplot(1,4,2);
imagesc(uint8(dimg_sm));
axis off;

ax = gca;
ax.Position = [0.255 0.00 0.24 1.00];

subplot(1,4,3);
imagesc(uint8(gt_sm));
axis off;

ax = gca;
ax.Position = [0.508 0.00 0.24 1.00];

subplot(1,4,4);
imagesc(uint8(dnn_img_sm));
axis off;

ax = gca;
ax.Position = [0.765 0.00 0.24 1.00];

print(plot_num, '-dpng', fullfile(save_path,strcat('dm_err_row1.png')));

plot_num = plot_num + 1;


figure(plot_num);
set(gcf,'position',([50,50,1300,300]),'color','w')

subplot(1,4,1);
image(gt_diff);
colormap(jet(max_depthmap_value+1));
axis off;

%max_error = max((gt_diff(:)));
max_error = max_depthmap_value;
cb = colorbar('fontweight','bold','FontSize', 12, 'Location', 'eastoutside');
cb.Label.String = 'Depth Map Error';
cb.Ticks = [0:5:max_error];
cb.TickLabels = num2str(cb.Ticks');
cb.Limits = [0 max_error];

ax = gca;
ax.Position = [0.00 0.00 0.24 1.00];
p = cb.Position;
p(2) = 0.03;
p(4) = 0.94;
cb.Position = p;

subplot('Position',[0.58 0.00 0.24 1.00]);
imagesc(uint8(gt_mask_over));
axis off;

subplot('Position',[0.325 0.00 0.24 1.00]);
imagesc(uint8(gt_mask_sm));
axis off;

print(plot_num, '-dpng', fullfile(save_path,strcat('dm_err_row2.png')));

plot_num = plot_num + 1;

return;

%% Plot the images and then save them

% plot the two rows
figure(plot_num); 
set(gcf,'position',([50,50,1200,250]),'color','w', 'Name','Combined');
imagesc(uint8(img_row1));
axis off;
plot_num = plot_num + 1;

figure(plot_num); 
set(gcf,'position',([50,50,1200,250]),'color','w', 'Name','Combined');
imagesc(uint8(img_row2));
axis off;
plot_num = plot_num + 1;

% save the two rows
if(save_data)
    imwrite(uint8(img_row1),fullfile(save_path,'dm_err_row1.png'));
    imwrite(uint8(img_row2),fullfile(save_path,'dm_err_row2.png'));
end
