format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% shape parameters

%color = {'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white', [0.5, 0.5, 0.5]};
% color = {[0 0 1]; [0 1 0]; [1 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; [0 0 0]; [1 1 1]; [0.5 0.5 0.5]};                    
% color_palette = 'basic';

% http://alumni.media.mit.edu/~wad/color/numbers.html
% color = {[0, 0, 0];[87, 87, 87]/255;[173, 35, 35]/255;[42, 75, 215]/255;...
%          [29, 105, 20]/255;[129, 74, 25]/255;[129, 38, 192]/255;[160, 160, 160]/255;...
%          [129, 197, 122]/255;[157, 175, 255]/255;[41, 208, 208]/255;[255, 146, 51]/255;...
%          [255, 238, 51]/255;[233, 222, 187]/255;[255, 205, 243]/255;[255, 255, 255]/255};
% color_palette = 'mit';

% 6-7-6 RGB color palette https://en.wikipedia.org/wiki/List_of_software_palettes
% green = [0, 42, 85, 128, 170, 212, 255];
% color = [];
% for r=0:5
%     for g=0:6
%         for b=0:5
%             color(end+1,:) = [51*r, green(g+1), 51*b]/255;
%         end
%     end
% end
% color_palette = '676';

% green 
% color = [35,44,41; 61,91,57; 113,114,80;  132,126,64]/255;
% color_palette = 'wood';

color = [];
color_palette = 'full';

commandwindow;

%% create the folders
save_path = 'D:/IUPUI/Test_data/tb19_test2/';

warning('off');
mkdir(save_path);
mkdir(save_path, 'images');
mkdir(save_path, 'depth_maps');
warning('on');

% the number of images to generate - not including the intensity variants
num_images = 100;
img_offset = 0;

%% load up the image generator

% lib_path = 'D:\Projects\simplex_noise\build\Release\';
% lib_name = 'sn_lib';
% hfile = 'D:\Projects\simplex_noise\include\sn_lib.h';
% 
% if(~libisloaded(lib_name))
%     [notfound, warnings] = loadlibrary(fullfile(lib_path, strcat(lib_name,'.dll')), hfile);
% end

 
% if(~libisloaded(lib_name))
%    fprintf('\nThe %s library did not load correctly!',  lib_name);    
% end


%% setup the blurring paramters

% gaussian kernel size
kernel_size = 61;

% sigma values for gaussian kernel
% this sigma ranges from 0 - 71 pixel blur radius:y = 0.1199x - 0.0984
sigma = [0.050,0.250,0.300,0.450,0.500,0.700,0.750,0.900,0.950,1.150,1.200,1.350,1.400,1.600,1.650,1.850,1.900,...
         2.100,2.150,2.300,2.350,2.550,2.600,2.800,2.850,3.050,3.100,3.300,3.350,3.500,3.550,3.750,3.800,4.000,...
         4.050,4.250,4.300,4.450,4.500,4.700,4.750,4.950,5.000,5.200,5.250,5.450,5.500,5.650,5.700,5.900,5.950,...
         6.150,6.200,6.400,6.450,6.650,6.700,6.850,6.900,7.100,7.150,7.350,7.400,7.650,7.700,7.900,7.950,8.250,...
         8.300,8.800,8.850];


% these are the candidate pixel blur radius to test
% 2 meters - 10 meters
% br1 = [9, 11, 13, 13, 14, 14, 15, 15, 15];
% br2 = [7, 4, 3, 2, 2, 1, 1, 1, 0];
% 30 - 49 meters - 1 meter increments
%br1 = [3, 4, 5, 6, 6, 7, 8, 9, 9, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15];
%br2 = [31, 30, 30, 29, 28, 27, 26, 26, 25, 24, 24, 23, 23, 22, 22, 21, 21, 20, 20, 19];
%tb16: 80-100
% br1 = [22, 23, 23, 24, 24, 25, 26, 26, 27, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31, 32];
% br2 = [45, 45, 44, 44, 43, 42, 42, 41, 41, 40, 40, 40, 39, 39, 38, 38, 37, 37, 37, 36, 36];
% tb17: 70-90
% br1 = [26,27,27,28,29,30,31,31,32,33,33,34,34,35,36,36,37,37,38,38,39];
% br2 = [41,40,39,38,38,37,36,35,35,34,33,33,32,31,31,30,30,29,29,28,28];
% tb18: 60-100
% br1 = [1,2,3,4,5,6,7,8,9,9,10,11,12,12,13,14,15,15,16,16,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,25,26,26,27];
% br2 = [26,25,24,23,22,21,21,20,19,18,17,16,16,15,14,14,13,12,12,11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,2,1,1];
% tb19: 60-100
br1 = [1,2,3,4,5,6,7,8,9,9,10,11,12,12,13,14,15,15,16,16,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,25,26,26,27];
br2 = [52,51,50,49,48,47,47,46,45,44,43,42,42,41,40,40,39,38,38,37,36,36,35,35,34,34,33,33,32,32,31,31,30,30,29,29,29,28,28,27,27];

% depth map values - arrange from lowest to highest with 0 being the lowest
% depthmap_range = [0:1:49];
depthmap_range = [0:1:(numel(br1)-1)];
max_depthmap = max(depthmap_range(:));
num_dm_values = numel(depthmap_range);

%% create all of the image generation parameters

% max number of depth map values for a single image
DM_N = floor(numel(depthmap_range)/4);

% initial number of objects at the first depthmap value
num_objects = 40;          

% block dimensions
blk_h = 50;
blk_w = 50;
max_blk_dim = max(blk_h, blk_w);

% build inmage dimensions
img_w = 400 + ceil(blk_w/2);
img_h = 400 + ceil(blk_h/2);
max_dim = max(img_w, img_h);

% make this cropping a mod 16 number
%img_w_range = 20:379;
%img_h_range = 20:379;
img_w_range = 17:400;
img_h_range = 17:400;

% intensity values to simulate different light conditions
%int_values = [0.2, 0.4, 0.6, 0.8, 1.0];
int_values = [1.0];

% x_min,x_max; y_min,y_max; min_r,max_r
rect = [ceil(max_blk_dim/7), ceil(max_blk_dim/5)];
circle = [ceil(max_blk_dim/7), ceil(max_blk_dim/5)];
polygon = [-ceil(max_blk_dim/5), ceil(max_blk_dim/5)];
shape_lims = {circle, polygon, rect};

rect_l = [ceil(max_dim/18), ceil(max_dim/14)];
circle_l= [ceil(max_dim/18), ceil(max_dim/14)];
polygon_l = [-ceil(max_dim/11), ceil(max_dim/11)];
shape_lims_l = {circle_l, polygon_l, rect_l};

% setting the image generation parameters for the simplex noise
% scale: 1 -> smallest, 0.0001 -> largest
octaves = 7;
persistence = 70/100;


%% start to create the images

fprintf('# %s\n\n', color_palette);
fprintf('%s\n\n', save_path);

% spmd
%     loadlibrary(fullfile(lib_path, strcat(lib_name,'.dll')), hfile);
%     
% end

tic;
parfor kdx=0:(num_images-1)
    
    % randomly generate DM_N depth values 
    if(depthmap_range(end) >= depthmap_range(1))
        D = randi([depthmap_range(1), depthmap_range(end)], 1, DM_N);
    else
        D = randi([depthmap_range(end), depthmap_range(1)], 1, DM_N);
    end
    D = sort(unique(D), 'descend');    
    
    % generate the first depthmap value - use the largest (farthest) value
%     dm = (max_depthmap/255)*ones(img_h, img_w, 3);
    dm = (D(1)/255)*ones(img_h, img_w, 3);
    
    for idx=2:numel(D)

        % get the number of shapes for a given depth map value
%         min_N = ceil(exp((3.7*(num_objects-idx))/num_objects));
%         max_N = ceil(exp((4.0*(num_objects-idx))/num_objects));
%         min_N = ceil(num_objects*(exp((4.0*((D(idx)-max_depthmap)))/max_depthmap)));
%         max_N = ceil(num_objects*(exp((1.0*((D(idx)-max_depthmap)))/max_depthmap)));

        min_N = ceil( ((num_objects)/(1+exp(-0.2*D(idx)+(0.1*num_objects))) ) + 3);
        max_N = ceil(1.25*min_N);
         
%         min_N = ceil(D(idx) + 5);
%         max_N = ceil(1.2*min_N);
        
        N = randi([min_N, max_N], 1);
        
        dm_blk = (D(idx)/255)*ones(blk_h, blk_w, 3);
        
        
        for jdx=1:N
            % the number of shapes in an image block
            S = randi([25,45], 1);
            %block = gen_rand_image_all(blk_h, blk_w, S, shape_lims);

            
            % generate random number to pick either the block or a circle
            shape_type = randi([0,1],1);
            switch(shape_type)
                case 0
                    A = randi([0, 89], 1, 1);
                    %block = imrotate(block, A, 'nearest', 'loose');
                    %blur_mask_blk_r = imrotate(blur_mask_blk, A, 'nearest', 'loose');

                    % check for a depthmap value of 0 and handle as a special case
                    if(D(idx) == 0)
                        dm_blk_r = imrotate(dm_blk+1, A, 'nearest', 'loose');
                        rotation_mask = dm_blk_r > 0;
                        dm_blk_r = imrotate(dm_blk, A, 'nearest', 'loose');                
                    else
                        dm_blk_r = imrotate(dm_blk, A, 'nearest', 'loose');
                        rotation_mask = dm_blk_r > 0;                
                    end
                case 1
                    R = floor(randi([floor(max_blk_dim*0.8), max_blk_dim],1)/2 + 0.5);
                    circle_mask = zeros(blk_h, blk_w);
                    circle_mask = insertShape(circle_mask, 'FilledCircle', [floor(blk_w/2), floor(blk_h/2), R], 'Color', [1,1,1], 'Opacity',1, 'SmoothEdges', false);
                    circle_mask = circle_mask(:,:,1);
                    %block = block.*circle_mask;
                    dm_blk_r = dm_blk.*circle_mask;
                    rotation_mask = circle_mask;
                    %blur_mask_blk_r = 1-circle_mask;
            end
            
            X = randi([1,img_w-blk_w], 1);
            Y = randi([1,img_h-blk_h], 1);
            
%             layer_img1 = overlay_with_mask(layer_img1, block, rotation_mask, X, Y);
%             layer_img2 = overlay_with_mask(layer_img2, block, rotation_mask, X, Y);
            dm = overlay_with_mask(dm, dm_blk_r, rotation_mask, X, Y);
%             blur_mask = overlay_with_mask(blur_mask, blur_mask_blk_r, rotation_mask, X, Y);
            
        end
        
%         % create the 2-D gaussian kernel using the depth map value indexing the sigma array
%         d_idx = D(idx) + 1;
% 
%         k1 = create_gauss_kernel(kernel_size, sigma( br1( d_idx ) + 1 ) );
%         k2 = create_gauss_kernel(kernel_size, sigma( br2( d_idx ) + 1 ) );
%         
%         % blur the layer and the blur_mask
%         LI_1 = imfilter(layer_img1, k1, 'corr', 'replicate', 'same');
%         BM_1 = imfilter(blur_mask, k1, 'corr', 'replicate', 'same');
%         LI_2 = imfilter(layer_img2, k2, 'corr', 'replicate', 'same');
%         BM_2 = imfilter(blur_mask, k2, 'corr', 'replicate', 'same');
%         
%         % bring the images back down to the original size
%         LI_1 = LI_1(1:img_h, 1:img_w, :);
%         BM_1 = BM_1(1:img_h, 1:img_w, :);
%         LI_2 = LI_2(1:img_h, 1:img_w, :);
%         BM_2 = BM_2(1:img_h, 1:img_w, :);
%         
%         % blending the current layer image and the previous image
%         img1 = (LI_1.*(1-BM_1) + (img1.*BM_1));
%         img2 = (LI_2.*(1-BM_2) + (img2.*BM_2));
        
        bp = 1;
        
    end


    dm = dm(img_h_range, img_w_range, :);
    
    image_num = num2str(kdx+img_offset, '%03d');
    dm_filename = strcat('depth_maps/dm_', image_num, '.png');
    imwrite(dm, strcat(save_path, dm_filename));
    
    for jdx=1:numel(int_values)
        
        % save the image file and depth maps
        image_int = num2str(int_values(jdx)*100, '%03d');
        
        img_filename1 = strcat('images/image_f1_', image_num, '.png');
        img_filename2 = strcat('images/image_f2_', image_num, '.png');
        
        fprintf('%s, %s, %s\n', img_filename1, img_filename2, dm_filename);
    end
    
end
toc;

%% save the image names to the standard image input file structure

save_name = strcat(save_path,'input_file_',datestr(now,'yyyymmdd_HHMMSS'),'.txt');
file_id = fopen(save_name, 'w');
fprintf(file_id, '# %s\n\n', color_palette);
fprintf(file_id, '%s\n\n', save_path);

for kdx=0:(num_images-1)
    % save the image file and depth maps
    image_num = num2str(kdx+img_offset, '%03d');
    dm_filename = strcat('depth_maps/dm_', image_num, '.png');
    
    for jdx=1:numel(int_values)
        image_int = num2str(int_values(jdx)*100, '%03d');
%         img_filename1 = strcat('images/image_f1_', image_num, '_', image_int, '.png');
%         img_filename2 = strcat('images/image_f2_', image_num, '_', image_int, '.png');
        img_filename1 = strcat('images/image_f1_', image_num, '.png');
        img_filename2 = strcat('images/image_f2_', image_num, '.png');
        fprintf(file_id, '%s, %s, %s\n', img_filename1, img_filename2, dm_filename);
    end
    
end

fclose(file_id);

