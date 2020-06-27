format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% create a surface map of the temp vs blur

startpath = 'D:\IUPUI\Test_Data\lens_raw\temp_test';
file_filter = {'*.png','PNG Files';'*.*','All Files' };

[image_files, image_path] = uigetfile(file_filter, 'Select Images to Analyze', startpath, 'MultiSelect', 'on');
if(image_path == 0)
    return;
end

save_path = 'D:\IUPUI\PhD\Images\Camera';

commandwindow;

%% read in the images

dct_thresh = 0;

img_row = 15;
num_img = numel(image_files);

cam_temp = cell(num_img,1);
img_fft = zeros(num_img,1); 
img_dct = zeros(num_img,1);
noise_fft = zeros(num_img,1);
noise_dct = zeros(num_img,1);

img_line = [];

center = [32,32];
radius = 11;

if(iscell(image_files))
    for idx=1:num_img
        
        [~, img_name,img_ext] = fileparts(image_files{idx}); 
        
        r_line = regexp(img_name,'_','split');
        cam_temp{idx,1} = r_line{4};
        cam_step{idx,1} = r_line{2};
        
        tmp_img = rgb2gray(imread(fullfile(image_path,image_files{idx})));
        
        mask = create_ring_mask(size(tmp_img), center, radius, 0);
        mask(center(2),center(1)) = 1;    % do this becasue of the math for removing the center radius
        
        mask_sum = sum(mask(:));
        
        img_noise = double(tmp_img).*(~mask);
        
        st = sum(img_noise(:));
        n = numel(img_noise) - sum(mask(:));
        st = floor(st/n);
        
        %img_noise = img_noise + (mask*(st/(numel(img_noise)-mask_sum)));
        img_noise = img_noise + (mask*st);
        
        %img_noise = tmp_img(:,[1:5,27:end]);
        
        img_fft(idx,1) = sum(sum(abs(fft2(tmp_img)/numel(tmp_img))));
        img_dct(idx,1) = get_dct_sharpness(tmp_img, dct_thresh);
        
        noise_fft(idx,1) = sum(sum(abs(fft2(img_noise)/numel(img_noise))));
        noise_dct(idx,1) = get_dct_sharpness(img_noise, dct_thresh);
        
        img{idx,1} = tmp_img;
        img_line(idx,:) = tmp_img(img_row,:);

    end
else
    [~, img_name,img_ext] = fileparts(image_files);
    img = imread(fullfile(image_path,image_files));
    img = img(crop_y1:crop_y2,crop_x1:crop_x2,:);
    
end


%% plot the lens over time

img_r = 3;
img_c = 5;

h = 10;
v = 10;

index = 1;
row = {};
for idx=1:img_r
    
    row{idx} = cat(2, img{index}, 255*ones(size(img{1},1),h), img{index+1}, 255*ones(size(img{1},1),h), ...
                      img{index+2}, 255*ones(size(img{1},1),h), img{index+3}, 255*ones(size(img{1},1),h), img{index+4});
                  
    imwrite(row{idx}, fullfile(save_path, strcat('lens_temp_row',num2str(idx),'.png')));

    index =  index + img_c;
end              
                      
combined_img = row{1};
for idx=2:img_r
    
    combined_img = cat(1, combined_img, 255*ones(v, size(row{1},2)), row{idx});
    
end

%imwrite(combined_img, fullfile(save_path, strcat('lens_temp_full.png')));


%% plot the sum of the fft

figure(plot_num)
set(gcf,'position',([50,50,1200,600]),'color','w')
hold on
box on
grid on

plot(1:num_img, img_fft, 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','b', 'MarkerSize',8);
plot(1:num_img, noise_fft, 'Color','r', 'LineStyle','-', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','r', 'MarkerSize',8);

set(gca,'FontSize', 13, 'fontweight','bold');

% X-Axis
xlim([1 num_img]);
xticks([1:1:num_img]);
%xtickformat('%2.1f');
%xticklabels(cam_step);
xticklabels(cam_temp);
xtickangle(90);
xlabel('Camera Temperature (^{o}C)', 'fontweight','bold','FontSize', 13);
%xlabel('Voltage Step', 'fontweight','bold','FontSize', 13);

% Y-Axis
%ylim([0, ceil(CoC_max/px_size)*px_size]);
ylim([0, ceil(max(img_fft)*.01)*100]);
%yticks(y_axis_ticks);
%yticklabels(y_axis_labels);
ylabel('Sum FFT Magnitude', 'fontweight','bold','FontSize', 13);

title(strcat('Image Sharpness vs. Camera Temperature'), 'fontweight','bold', 'FontSize',16);
%title(strcat('Lens Sharpness vs. Voltage Step'), 'fontweight','bold', 'FontSize',16);

lgd = legend('Image Energy','Noise Energy', 'location','northwest', 'orientation','vertical');

ax = gca;
ax.Position = [0.07 0.16 0.91 0.77];

print(plot_num, '-dpng', fullfile(save_path,strcat('lens_temp_fft_v2.png')));

plot_num = plot_num + 1;

%% plot the dct values

noise_diff = (img_dct-noise_dct);

figure(plot_num)
%set(gcf,'position',([50,50,1200,600]),'color','w')
set(gcf,'position',([50,50,700,500]),'color','w')
hold on
box on
grid on

plot(1:num_img, img_dct, 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','o', 'MarkerFaceColor','b', 'MarkerSize',4);
plot(1:num_img, noise_dct, 'Color','r', 'LineStyle','-', 'LineWidth',1, 'Marker','d', 'MarkerFaceColor','r', 'MarkerSize',4);
%plot(1:num_img, noise_diff, 'Color','g', 'LineStyle','-.', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','g', 'MarkerSize',8);

set(gca,'FontSize', 13, 'fontweight','bold');

% X-Axis
xlim([1 num_img]);
xticks([1:1:num_img]);
%xtickformat('%2.1f');
%xticklabels(cam_step);
xticklabels(cam_temp);
xtickangle(90);
xlabel('Camera Temperature (^{o}C)', 'fontweight','bold','FontSize', 13);
%xlabel('Voltage Step', 'fontweight','bold','FontSize', 13);

% Y-Axis
% y_min = floor(min(cat(1,img_dct(:), noise_dct(:), noise_diff(:)))*10)/10;
% y_max = ceil(max(cat(1,img_dct(:), noise_dct(:), noise_diff(:)))*10)/10;

y_min = floor(min(cat(1,img_dct(:), noise_dct(:)))*10)/10;
y_max = ceil(max(cat(1,img_dct(:), noise_dct(:)))*10)/10;

ylim([y_min, y_max]);
%yticks(y_axis_ticks);
%yticklabels(y_axis_labels);
ylabel('DCT Sharpness Measure', 'fontweight','bold','FontSize', 13);

title(strcat('Image Sharpness vs. Camera Temperature'), 'fontweight','bold', 'FontSize',16);
%title(strcat('Lens Sharpness vs. Voltage Step'), 'fontweight','bold', 'FontSize',16);

lgd = legend({'Image', 'Noise'}, 'location','northwest', 'orientation','vertical');

ax = gca;
%ax.Position = [0.07 0.16 0.91 0.77];        % 1200,600
ax.Position = [0.09 0.165 0.88 0.77];      % 700, 500

print(plot_num, '-dpng', fullfile(save_path,strcat('lens_temp_dct_v3.png')));

plot_num = plot_num + 1;
