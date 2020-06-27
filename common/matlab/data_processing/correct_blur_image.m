format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% get the path to save the results
startpath = 'D:\IUPUI\Test_Data\';
save_path = uigetdir(startpath,'Select Save Folder');

if(save_path == 0)
    return;
end

%% load in the mat file
% load the mat file with the kernel used to blur the image
% - kernel: contains the kernels used to blur the defocused image
% - ring mask params: start, step, end, center

% need to offset the center point for the images
% Original images are 728 x 736 (hxw)
% but the converted images are 630x630 (hxw)
% centered -> [49 - 630 - 49] x [53 - 630 - 53] (hxw)
% ring center = 306,196 (x,y) -> 253, 147 (x,y)

file_filter = {'*.mat','MAT-files (*.mat)'; '*.*','All Files' };

% select the mat file that contains the kernel and mask info
startpath = 'D:\IUPUI\Test_Data\rw_raw4\k64';
[kernel_file, data_path] = uigetfile(file_filter, 'Select MATLAB Data Kernel File', startpath);

if(data_path == 0)
    return;
end

% this file shold contain the following:
% - kernel
% - ring_center = [306,178];    %306,196 without the 18 removed from the top
% - r_step = 20;
% - r_outer = 80:r_step:max(img_w,img_h);

load(fullfile(data_path, kernel_file));

%% get the config file for the raw data

file_filter = {'*.txt','Text Files';'*.*','All Files' };
startpath = 'D:\IUPUI\Test_Data\'; 

[cfg_file, cfg_file_path] = uigetfile(file_filter, 'Select Configuration File', startpath, 'MultiSelect', 'on');
if(cfg_file_path == 0)
    return;
end

v_step = [128];
commandwindow;

%% 
if(iscell(cfg_file))
    for idx=1:numel(cfg_file)
        process_blur_correction(cfg_file{idx}, cfg_file_path(1:end-1), v_step, save_path, kernel, mask_params);
        fprintf('---------------------------------------------------------\n');
    end
else
    process_blur_correction(cfg_file, cfg_file_path(1:end-1), v_step, save_path, kernel, mask_params);
end

% cfg_parts = strsplit(cfg_file,{'_','.'});
% 
% scenario_name = cfg_parts{2};
% camera_side = cfg_parts{3};
% 
% folders = strsplit(cfg_file_path,filesep);
% data_path = fullfile(folders{1:end-2},scenario_name);
% 
% %% parse the contents of the file
% 
% [match_params] = parse_input_parameters(fullfile(cfg_file_path, cfg_file));
% 
% img_x_off = str2double(match_params{1}{1});
% img_y_off = str2double(match_params{1}{2});
% img_w = str2double(match_params{1}{3});
% img_h = str2double(match_params{1}{4});
% img_rot = str2double(match_params{1}{5});
% 
% img_crop_w = [img_x_off:img_x_off+img_w-1];
% img_crop_h = [img_y_off:img_y_off+img_h-1];
% 
% % get the exposre directory listing
% exp_listing = dir(strcat(data_path,filesep,camera_side,filesep,'exp*'));
%     
% for idx=1:numel(exp_listing)
%     
%     exp_name = exp_listing(idx).name;
%     img_listing=[];
%     img_listing = dir(strcat(exp_listing(idx).folder, filesep, exp_name, filesep,'*',num2str(v_step(1)),'*.png'));
%     
%     img_parts = strsplit(img_listing.name, {'_'});
%     
%     defocus_img = double(imread(fullfile(img_listing(1).folder,img_listing(1).name)));
%     
%     [corr_img] = apply_blur_correction(defocus_img, kernel, mask_params);
%     
%     save_img_name = strcat(img_parts{1}, '_', img_parts{2}, '_', img_parts{3}, '_corr.png');
%     
%     if(img_rot ~= 0.0)
%         img_r = imrotate(corr_img, img_rot, 'bilinear');
%         img_cr = img_r(img_crop_h, img_crop_w, :);
%     else
%         img_cr = corr_img(img_crop_h, img_crop_w, :);
%     end
% 
%     save_location = fullfile(save_path, scenario_name, camera_side, exp_name);
%     warning('off','all')
%     mkdir(save_location);
%     fprintf('Saving: %s\n', fullfile(save_location, save_img_name));
%     imwrite(uint8(img_cr), fullfile(save_location, save_img_name));
%     warning('on','all')
% end

fprintf('Complete!\n');

