% batch crop images

format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% get the images
startpath = 'D:\IUPUI\Test_Data\lens_raw\temp_test';
file_filter = {'*.png','PNG Files';'*.*','All Files' };

[image_files, image_path] = uigetfile(file_filter, 'Select Files to Crop', startpath, 'MultiSelect', 'on');
if(image_path == 0)
    return;
end

%save_path = 'D:\IUPUI\Test_Data\lens_raw\temp_test\crop';
%save_path = 'D:\IUPUI\Test_Data\lens_raw\nrzu\crop';
% save_path = 'D:\IUPUI\Test_Data\lens_raw\nrzd\crop';
save_path = image_path;

commandwindow;

%%


crop_x1 = 1;
crop_y1 = 1;
crop_x2 = 384;
crop_y2 = 384;

% temp crop
% crop_x1 = 363;  %379;
% crop_y1 = 286;  %302;
% crop_x2 = crop_x1+63;
% crop_y2 = crop_y1+63;

% crop_x1 = 20;
% crop_y1 = 60;
% crop_x2 = crop_x1+599;
% crop_y2 = crop_y1+599;

% nrz tests
% crop_x1 = 38;
% crop_y1 = 35;
% crop_x2 = crop_x1+127;
% crop_y2 = crop_y1+127;


figure(plot_num);
if(iscell(image_files))
    for idx=1:numel(image_files)
        
        [~, img_name,img_ext] = fileparts(image_files{idx});        
        img = imread(fullfile(image_path,image_files{idx}));
        img = img(crop_y1:crop_y2,crop_x1:crop_x2,:);
        
        image(img);
        fprintf('Writing: %s\n',fullfile(image_path,strcat(img_name,'_crop.png')));
        imwrite(img,fullfile(save_path,strcat(img_name,'_crop.png')));
        pause(0.5);
    end
else
    [~, img_name,img_ext] = fileparts(image_files);
    img = imread(fullfile(image_path,image_files));
    img = img(crop_y1:crop_y2,crop_x1:crop_x2,:);

    image(img);
    fprintf('Writing: %s\n',fullfile(image_path,strcat(img_name,'_crop.png')));
    imwrite(img,fullfile(save_path,strcat(img_name,'_crop.png')));
    
end

%%
fprintf('Complete!\n');
