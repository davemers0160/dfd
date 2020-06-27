format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;


%% get the directory for the images

img_path = uigetdir(startpath, 'Select Folder with Images');
image_ext = '*.png';

if(img_path == 0)
    return;
end

listing = dir(strcat(img_path, '/', image_ext));

save_dir = img_path;

%% run through each image in the list
commandwindow;

% this defines the expected maximum blur radius
max_blur_radius = 64;

% this is a value toprovide a small buffer to exceed before being counted
offset = 1/255;
    
fprintf('File Name\t\t\t\t# of Pixels Blurred\n')
fprintf('-----------------------------------------------------\n')
for idx=1:numel(listing)

    % load in an image and get its size
    img_file = fullfile(listing(idx).folder, listing(idx).name);
    img = double(imread(img_file));
    [img_h, img_w, img_c] = size(img);

    % find the rough center points assuming that the knife edge is right to
    % left and that the darker portion is on the left
    img_cw = floor(img_w/2);    
    width_range = max(0,img_cw-max_blur_radius):1:min(img_w,img_cw+max_blur_radius);
    
    % just use a single line to determine the blur amount
    img_line = img(floor(img_h/2), width_range, 1);

    % find the areas where limits are met
    low_limit = min(img_line(:));
    high_limit = max(img_line(:));
    
    match = (img_line > (low_limit+offset)) == (img_line < (high_limit-offset));
    num = sum(match);
    
    
    fprintf('%03d: %s, \t%02d\n', (idx-1), listing(idx).name, num);

%     figure(1)
%     plot(img_line, 'b');
%     hold on;
%     plot(match*max(img_line(:)), 'r');
%     hold off;

end

fprintf('-----------------------------------------------------\n')



