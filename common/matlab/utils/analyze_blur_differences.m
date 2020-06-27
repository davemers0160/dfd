format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[scriptpath,  filename, ext] = fileparts(full_path);
%mkdir(strcat(scriptpath,'\Images'));
plot_num = 1;

%% Open the images

save_path = 'D:\IUPUI\PhD\Images\datasets';

startpath = 'D:\IUPUI\Test_Data\rw4\';
file_filter = {'*.png','PNG Files';'*.*','All Files' };

[infocus_file, infocus_path] = uigetfile(file_filter, 'Select In-focus Image to Analyze', startpath, 'MultiSelect', 'off');

if(infocus_path == 0)
    return;
end


[defocus_file, defocus_path] = uigetfile(file_filter, 'Select Out-of-Focus Image to Analyze', infocus_path, 'MultiSelect', 'off');

if(defocus_path == 0)
    return;
end



% if((numel(image_path)==1 && (image_path == 0)) || ~iscell(image_files))
%     return;
% end

commandwindow;

fprintf('In-focus file:     %s\n', infocus_file);
fprintf('Out-of-focus file: %s\n', defocus_file);


%% get the image info

folders = regexp(infocus_path, filesep, 'split');
scene{1} = folders{end-3};

folders = regexp(defocus_path, filesep, 'split');
scene{2} = folders{end-3};

r = 250;
name = {'In-Focus', 'Out-of-Focus'};


img{1} = double(rgb2gray(imread(fullfile(infocus_path, infocus_file))));
img{2} = double(rgb2gray(imread(fullfile(defocus_path, defocus_file))));


for idx=1:numel(img)
    
    figure(plot_num);
    imagesc(2*abs(img{idx}));
    colormap(gray(256));
    title(strcat(scene{1}, 32, name{idx}));
    plot_num = plot_num + 1;
end

%% plot 

figure(plot_num);
hold on; 
grid on;
box on;
plot(img{1}(r,:), 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','b', 'MarkerSize',8); 
plot(img{2}(r,:), 'Color','r', 'LineStyle','--', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','b', 'MarkerSize',8);
title(strcat(scene{1}, 32, '-', 32, num2str(r)));
plot_num = plot_num + 1;


%% do the FFT of the image

img_fft = {};
for idx=1:numel(img)
    img_fft{idx} = fftshift(fft2(img{idx})/numel(img{idx}(:)));
    figure(plot_num);
    imagesc(2*abs(img_fft{idx}),[0,20]);
    colormap(jet(256));
    title(strcat('FFT:', 32, scene{1}, 32, name{idx}));
    plot_num = plot_num + 1;
end

%% plot the FFT

figure(plot_num);
hold on; 
grid on;
box on;
plot(2*abs((img_fft{1}(r,:))), 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','b', 'MarkerSize',8); 
plot(2*abs((img_fft{2}(r,:))), 'Color','r', 'LineStyle','--', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','b', 'MarkerSize',8);
title(strcat('FFT:', 32, scene{1}, 32, '-', 32, num2str(r)));
plot_num = plot_num + 1;




