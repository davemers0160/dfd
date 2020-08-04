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
max_blur_radius = 150;

sk = create_1D_gauss_kernel(5, 1.0);

mf = [1 1 1];

dx = [-0.5 0 0.5];

fx = [1 -2 1];

% this is a value toprovide a small buffer to exceed before being counted
offset = 1/255;
    
fprintf('File Name\t\t\t\t# of Pixels Blurred\n')
fprintf('-----------------------------------------------------\n')
for idx=1:numel(listing)

    % load in an image and get its size
    img_file = fullfile(listing(idx).folder, listing(idx).name);
    img = double(rgb2gray(imread(img_file)));
    [img_h, img_w, img_c] = size(img);

    % find the rough center points assuming that the knife edge is right to left
    img_s = img(floor(img_h/2-5:img_h/2+5),:);
    img_line = mean(img_s, 1);
    img_cw = floor(img_w/2);
    width_range = max(0,img_cw-max_blur_radius):1:min(img_w,img_cw+max_blur_radius);
    
    % just use a single line to determine the blur amount
    %img_line = conv(img_line, sk, 'same');
    img_line = img_line(width_range);
    
    % find the areas where limits are met
    low_limit = ceil(min(img_line(:))) + 1;
    high_limit = floor(max(img_line(:))) - 1;
    mid_limit = (high_limit + low_limit)/2;
    
    % find the point where the slope is changing
    % we assume that there enough distance between the min and max along
    % with noise
    %index = find(
    [~, mid_idx] = min(abs(img_line - mid_limit));
    
    
    % get the laplacian to see where the inflection point is located
    lap = conv(img_line, fx,'same');
    lap = lap(2:end-1);
    [min_lap, min_idx] = min(lap);
    [max_lap, max_idx] = max(lap);
    
    dx1 = diff(img_line, 1);
    dx2 = diff(img_line, 2);

    
    match = (img_line > (low_limit+offset)) == (img_line < (high_limit-offset));
    %match = conv(match, mf, 'same');
    match = bwareafilt(match,1);

    num = sum(match);
    
    fprintf('%03d: %s, \t%02d\n', (idx-1), listing(idx).name, num);

    figure(1)
    plot(img_line, 'b');
    hold on;
    plot(match*max(img_line(:)), 'r');
    %plot(img_line2, 'g');
    hold off;

end

fprintf('-----------------------------------------------------\n')



