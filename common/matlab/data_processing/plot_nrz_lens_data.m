
format long g
format compact
clc
close all
clearvars

plot_num = 1;

%% load the data

startpath = 'D:\IUPUI\Test_Data\lens_raw\nrzu';
nrzu_data_path = uigetdir(startpath, 'Select NRZU Crop Folder');
if(nrzu_data_path == 0)
    return;
end

startpath = 'D:\IUPUI\Test_Data\lens_raw\nrzd';
nrzd_data_path = uigetdir(startpath, 'Select NRZD Crop Folder');
if(nrzd_data_path == 0)
    return;
end

save_path = 'D:\IUPUI\PhD\Images\Camera';

commandwindow;

% nrzu = load('D:\IUPUI\Test_Data\lens_raw\nrzu\nrzu_data.mat');
% nrzd = load('D:\IUPUI\Test_Data\lens_raw\nrzd\nrzd_data.mat');

%% cycle through the folders 
center = [65,62];
radius = 12;
dct_thresh = 0.0;

nrzu_listing = dir(strcat(nrzu_data_path,filesep,'*.png'));
num_nrzu_trials = numel(nrzu_listing);

cam_step = cell(num_nrzu_trials,1);
nrzu_img_fft = zeros(num_nrzu_trials,1); 
nrzu_noise_fft = zeros(num_nrzu_trials,1);
nrzu_img_dct = zeros(num_nrzu_trials,1); 
nrzu_noise_dct = zeros(num_nrzu_trials,1);

for idx=1:num_nrzu_trials
    
        [~, img_name,img_ext] = fileparts(nrzu_listing(idx).name); 
        
        r_line = regexp(img_name,'_','split');
        cam_step{idx,1} = r_line{2};
        
        tmp_img = rgb2gray(imread(fullfile(nrzu_listing(idx).folder,nrzu_listing(idx).name)));
        
        mask = create_ring_mask(size(tmp_img), center, radius, 0);
        mask(center(2),center(1)) = 1;    % do this becasue of the math for removing the center radius
        
        mask_sum = sum(mask(:));
        
        nrzu_noise = double(tmp_img).*(~mask);
        
        st = sum(nrzu_noise(:));
        n = numel(nrzu_noise) - sum(mask(:));
        st = floor(st/n);
        
        nrzu_noise = nrzu_noise + (mask*st);
        
        %img_noise = tmp_img(:,[1:5,27:end]);
        
        nrzu_img_fft(idx,1) = sum(sum(abs(fft2(tmp_img)/numel(tmp_img))));
        nrzu_noise_fft(idx,1) = sum(sum(abs(fft2(nrzu_noise)/numel(nrzu_noise))));  
        
        nrzu_img_dct(idx,1) = get_dct_sharpness(tmp_img, dct_thresh);
        nrzu_noise_dct(idx,1) = get_dct_sharpness(nrzu_noise, dct_thresh);
        
end

nrzd_listing = dir(strcat(nrzd_data_path,filesep,'*.png'));
num_nrzd_trials = numel(nrzd_listing);

nrzd_img_fft = zeros(num_nrzd_trials,1); 
nrzd_noise_fft = zeros(num_nrzd_trials,1);
nrzd_img_dct = zeros(num_nrzd_trials,1); 
nrzd_noise_dct = zeros(num_nrzd_trials,1);

for idx=1:num_nrzd_trials
    
        tmp_img = rgb2gray(imread(fullfile(nrzd_listing(idx).folder,nrzd_listing(idx).name)));
        
        mask = create_ring_mask(size(tmp_img), center, radius, 0);
        mask(center(2),center(1)) = 1;    % do this becasue of the math for removing the center radius
        
        mask_sum = sum(mask(:));
        
        nrzd_noise = double(tmp_img).*(~mask);
        
        st = sum(nrzd_noise(:));
        n = numel(nrzd_noise) - sum(mask(:));
        st = floor(st/n);
        
        nrzd_noise = nrzd_noise + (mask*st);
        
        %img_noise = tmp_img(:,[1:5,27:end]);
        
        nrzd_img_fft(idx,1) = sum(sum(abs(fft2(tmp_img)/numel(tmp_img))));
        nrzd_noise_fft(idx,1) = sum(sum(abs(fft2(nrzd_noise)/numel(nrzd_noise))));  
        
        nrzd_img_dct(idx,1) = get_dct_sharpness(tmp_img, dct_thresh);
        nrzd_noise_dct(idx,1) = get_dct_sharpness(nrzd_noise, dct_thresh);
end


%% plot the data fft version

figure(plot_num)
set(gcf,'position',([50,50,1200,600]),'color','w')
hold on
box on
grid on

% plot the data
plot(1:num_nrzu_trials, nrzu_img_fft, 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','o', 'MarkerFaceColor','b', 'MarkerSize',4);
plot(1:num_nrzd_trials, nrzd_img_fft, 'Color','r', 'LineStyle','-', 'LineWidth',1, 'Marker','d', 'MarkerFaceColor','r', 'MarkerSize',4);

%plot the noise
% plot(1:num_nrzu_trials, nrzu_noise_fft, 'Color','b', 'LineStyle','--', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','b', 'MarkerSize',8);
% plot(1:num_nrzd_trials, nrzd_noise_fft, 'Color','r', 'LineStyle','--', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','r', 'MarkerSize',8);

set(gca,'FontSize', 13, 'fontweight','bold');

% X-Axis
xlim([1 num_nrzu_trials]);
xticks([1:1:num_nrzu_trials]);
%xtickformat('%2.1f');
xticklabels(cam_step);
xtickangle(90);
xlabel('Voltage Step', 'fontweight','bold','FontSize', 13);

% Y-Axis
%ylim([0, ceil(CoC_max/px_size)*px_size]);
ylim([0, ceil(max(nrzu_img_fft)*.01)*100]);
%yticks(y_axis_ticks);
%yticklabels(y_axis_labels);
ylabel('Sum FFT Magnitude', 'fontweight','bold','FontSize', 13);

title(strcat('Lens Sharpness vs. Voltage Step'), 'fontweight','bold', 'FontSize',16);

%legend('127 \rightarrow 143','143 \rightarrow 127', '127 \rightarrow 143 Noise','143 \rightarrow 127 Noise','location','northwest', 'orientation','vertical')
legend({'127 \rightarrow 143','143 \rightarrow 127'}, 'location','northwest', 'orientation','vertical')

ax = gca;
ax.Position = [0.07 0.11 0.91 0.82];

print(plot_num, '-dpng', fullfile(save_path,strcat('nrz_lens_step_fft.png')));

plot_num = plot_num + 1;

%% plot the dct version

figure(plot_num)
%set(gcf,'position',([50,50,1200,600]),'color','w')
set(gcf,'position',([50,50,700,500]),'color','w')
hold on
box on
grid on

% plot the data
plot(1:num_nrzu_trials, nrzu_img_dct, 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','o', 'MarkerFaceColor','b', 'MarkerSize',4);
plot(1:num_nrzd_trials, nrzd_img_dct, 'Color','r', 'LineStyle','--', 'LineWidth',1, 'Marker','d', 'MarkerFaceColor','r', 'MarkerSize',4);

%plot the noise
% plot(1:num_nrzu_trials, nrzu_noise_dct, 'Color','b', 'LineStyle','--', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','b', 'MarkerSize',8);
% plot(1:num_nrzd_trials, nrzd_noise_dct, 'Color','r', 'LineStyle','--', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','r', 'MarkerSize',8);

set(gca,'FontSize', 13, 'fontweight','bold');

% X-Axis
xlim([1 num_nrzu_trials]);
xticks([1:1:num_nrzu_trials]);
%xtickformat('%2.1f');
xticklabels(cam_step);
xtickangle(90);
xlabel('Voltage Step', 'fontweight','bold','FontSize', 13);

% Y-Axis
% y_max = ceil(max(max(cat(1,nrzu_img_dct(:),nrzd_img_dct(:))), max(cat(1,nrzu_noise_dct(:),nrzd_noise_dct(:))))*10)/10;
% y_min = floor(min(min(cat(1,nrzu_img_dct(:),nrzd_img_dct(:))), min(cat(1,nrzu_noise_dct(:),nrzd_noise_dct(:))))*10)/10;
y_max = ceil(max(max(cat(1,nrzu_img_dct(:),nrzd_img_dct(:))))*10)/10;
y_min = floor(min(min(cat(1,nrzu_img_dct(:),nrzd_img_dct(:))))*20)/20;

ylim([y_min, y_max]);
%yticks(y_axis_ticks);
%yticklabels(y_axis_labels);
ytickformat('%1.2f');
ylabel('DCT Sharpness Measure', 'fontweight','bold','FontSize', 13);

title(strcat('Image Sharpness vs. Voltage Step'), 'fontweight','bold', 'FontSize',16);

%legend('127 \rightarrow 143','143 \rightarrow 127', '127 \rightarrow 143 Noise','143 \rightarrow 127 Noise','location','northwest', 'orientation','vertical')
legend({'127 \rightarrow 143','143 \rightarrow 127'}, 'location','northwest', 'orientation','vertical')

ax = gca;
%ax.Position = [0.07 0.11 0.91 0.82];       % 1200, 600
ax.Position = [0.105 0.125 0.87 0.81];      % 700, 500

print(plot_num, '-dpng', fullfile(save_path,strcat('nrz_lens_step_dct_v2.png')));

plot_num = plot_num + 1;

