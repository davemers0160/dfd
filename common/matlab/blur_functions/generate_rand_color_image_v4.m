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
save_path = 'D:/IUPUI/Test_data/tb22d_test/';

warning('off');
mkdir(save_path);
mkdir(save_path, 'images');
mkdir(save_path, 'depth_maps');
warning('on');

% the number of images to generate - not including the intensity variants
num_images = 30;
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
kernel_size = 79;

% sigma values for gaussian kernel
% this sigma ranges from 0 - 71 pixel blur radius:y = 0.1199x - 0.0984
% sigma = [0.050,0.250,0.300,0.450,0.500,0.700,0.750,0.900,0.950,1.150,1.200,1.350,1.400,1.600,1.650,1.850,1.900,...
%          2.100,2.150,2.300,2.350,2.550,2.600,2.800,2.850,3.050,3.100,3.300,3.350,3.500,3.550,3.750,3.800,4.000,...
%          4.050,4.250,4.300,4.450,4.500,4.700,4.750,4.950,5.000,5.200,5.250,5.450,5.500,5.650,5.700,5.900,5.950,...
%          6.150,6.200,6.400,6.450,6.650,6.700,6.850,6.900,7.100,7.150,7.350,7.400,7.650,7.700,7.900,7.950,8.250,...
%          8.300,8.800,8.850];

% new sigma with fixes for the double to uint8 conversion to images
% sigma range 0 - 64 pixel blur radius: y = 0.17311101398602x - 0.0498139860140584
% sigma = [0.001,0.284,0.285,0.586,0.587,0.912,0.913,1.248,1.249,1.588,1.589,1.930,1.931,2.274,2.275,2.619,2.620,...
%          2.964,2.965,3.309,3.310,3.655,3.656,4.001,4.002,4.347,4.348,4.694,4.695,5.040,5.041,5.387,5.388,5.733,...
%          5.734,6.080,6.081,6.426,6.427,6.773,6.774,7.120,7.121,7.466,7.467,7.813,7.814,8.160,8.161,8.507,8.508,...
%          8.855,8.856,9.205,9.206,9.556,9.557,9.912,9.913,10.274,10.275,10.647,10.648,11.035,11.036];
sigma = [0.0005,0.1725,0.3450,0.5175,0.6900,0.8625,1.0350,1.2075,1.3800,1.5525,1.7250,1.8975,2.0700,2.2425,2.4150,...
         2.5875,2.7600,2.9325,3.1050,3.2775,3.4500,3.6225,3.7950,3.9675,4.1400,4.3125,4.4850,4.6575,4.8300,5.0025,...
         5.1750,5.3475,5.5200,5.6925,5.8650,6.0375,6.2100,6.3825,6.5550,6.7275,6.9000,7.0725,7.2450,7.4175,7.5900,...
         7.7625,7.9350,8.1075,8.2800,8.4525,8.6250,8.7975,8.9700,9.1425,9.3150,9.4875,9.6600,9.8325,10.0050,10.1775,....
         10.3500,10.5225,10.6950,10.8675,11.0400,11.2125,11.3850,11.5575,11.7300,11.9025,12.0750,12.2475];


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
% tb19/21/22: 60-100
% br1 = [1,2,3,4,5,6,7,8,9,9,10,11,12,12,13,14,15,15,16,16,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,25,26,26,27];
% br2 = [52,51,50,49,48,47,47,46,45,44,43,42,42,41,40,40,39,38,38,37,36,36,35,35,34,34,33,33,32,32,31,31,30,30,29,29,29,28,28,27,27];
% tb20/21a/22a: 60-80
% br1 = [1,2,3,4,5,6,7,8,9,9,10,11,12,12,13,14,15,15,16,16,17];
% br2 = [17,16,15,14,13,12,12,11,10,9,8,8,7,6,5,5,4,3,3,2,1];
% tb20/21b/22b: 80-100
% br1 = [1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9,10,10,11];
% br2 = [11,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,3,2,2,1,1];
% tb22c: 60-100 V
br1 = [1,2,3,4,5,6,7,8,9,9,10,11,12,12,13,14,15,15,16,16,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,25,26,26,27];
br2 = [26,25,24,23,22,21,21,20,19,18,17,16,16,15,14,14,13,12,12,11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,2,1,1];

% depth map values - arrange from lowest to highest with 0 being the lowest
% depthmap_range = [0:1:49];
depthmap_range = [0:1:(numel(br1)-1)];
max_depthmap = max(depthmap_range(:));
num_dm_values = numel(depthmap_range);

%% create all of the image generation parameters

% max number of depth map values for a single image
DM_N = floor(numel(depthmap_range)/2);

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
for kdx=0:(num_images-1)
    
    img1 = [];
    % create an image as a background instead of a solid color
    if(strcmp(color_palette, 'full'))
        img1 = gen_rand_image_all(img_h, img_w, 450, shape_lims_l);
    else
        seed = int32(double(intmin('int32')) + double(intmax('uint32'))*rand(1));
        calllib(lib_name, 'init', seed);
        scale = 1;
        for r=1:img_h
            for c=1:img_w
                index = calllib(lib_name, 'octave_evaluate', r, c, scale, octaves, persistence);
                img1(r,c,:) = color(index+1, :);
            end
        end
    end
    img2 = img1;
    
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
    
    % blur that backgraound according to the depthmap value
    k1 = create_gauss_kernel(kernel_size, sigma( br1( D(1) + 1 ) + 1 ) );
    k2 = create_gauss_kernel(kernel_size, sigma( br2( D(1) + 1 ) + 1 ) );

    % blur the layer and the blur_mask
    img1 = imfilter(img1, k1, 'corr', 'replicate', 'same');
    img2 = imfilter(img2, k2, 'corr', 'replicate', 'same');
   
    for idx=2:numel(D)
        
        layer_img1 = img1;
        layer_img2 = img2;        
        
        % get the number of shapes for a given depth map value
%         min_N = ceil(exp((3.7*(num_objects-idx))/num_objects));
%         max_N = ceil(exp((4.0*(num_objects-idx))/num_objects));
%         min_N = ceil(num_objects*(exp((4.0*((D(idx)-max_depthmap)))/max_depthmap)));
%         max_N = ceil(num_objects*(exp((1.0*((D(idx)-max_depthmap)))/max_depthmap)));
        %min_N = ceil( ((num_objects)/(1+exp(-0.2*D(idx)+(0.1*num_objects))) ) + 3);
        min_N = ceil( ((max_depthmap)/(1+exp(-0.2*D(idx)+(0.1*max_depthmap))) ) + 3);
        max_N = ceil(1.25*min_N);
        
        N = randi([min_N, max_N], 1);
        
        dm_blk = (D(idx)/255)*ones(blk_h, blk_w, 3);
        
        % create the overlay mask and mask block
        blur_mask = ones(img_h+blk_h, img_w+blk_w);
        blur_mask_blk = zeros(blk_h, blk_w);
        
        blur_mask_blk_r = [];
        dm_blk_r = [];
        rotation_mask = [];
        
        scale = 1/(max_depthmap + 1 - D(idx));
        
        for jdx=1:N
            if(strcmp(color_palette, 'full'))
                % the number of shapes in an image block
                S = randi([25,45], 1);
                block = gen_rand_image_all(blk_h, blk_w, S, shape_lims);
            else
                block = zeros(blk_h, blk_w, 3);            
                seed = int32(double(intmin('int32')) + double(intmax('uint32'))*rand(1));
                calllib(lib_name, 'init', seed);
                for r=1:blk_h
                    for c=1:blk_w
                        index = calllib(lib_name, 'octave_evaluate', r, c, scale, octaves, persistence);
                        block(r,c,:) = color(index+1, :);
                    end
                end
            end
            
            % generate random number to pick either the block or a circle
            shape_type = randi([0,1],1);
            switch(shape_type)
                case 0
                    A = randi([0, 89], 1, 1);
                    block = imrotate(block, A, 'nearest', 'loose');
                    blur_mask_blk_r = imrotate(blur_mask_blk, A, 'nearest', 'loose');

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
                    block = block.*circle_mask;
                    dm_blk_r = dm_blk.*circle_mask;
                    rotation_mask = circle_mask;
                    blur_mask_blk_r = 1-circle_mask;
            end
            
            X = randi([1,img_w-blk_w], 1);
            Y = randi([1,img_h-blk_h], 1);
            
            layer_img1 = overlay_with_mask(layer_img1, block, rotation_mask, X, Y);
            layer_img2 = overlay_with_mask(layer_img2, block, rotation_mask, X, Y);
            dm = overlay_with_mask(dm, dm_blk_r, rotation_mask, X, Y);
            blur_mask = overlay_with_mask(blur_mask, blur_mask_blk_r, rotation_mask, X, Y);
            
        end
        
        % create the 2-D gaussian kernel using the depth map value indexing the sigma array
        d_idx = D(idx) + 1;

        k1 = create_gauss_kernel(kernel_size, sigma( br1( d_idx ) + 1 ) );
        k2 = create_gauss_kernel(kernel_size, sigma( br2( d_idx ) + 1 ) );
        
        % blur the layer and the blur_mask
        LI_1 = imfilter(layer_img1, k1, 'corr', 'replicate', 'same');
        BM_1 = imfilter(blur_mask, k1, 'corr', 'replicate', 'same');
        LI_2 = imfilter(layer_img2, k2, 'corr', 'replicate', 'same');
        BM_2 = imfilter(blur_mask, k2, 'corr', 'replicate', 'same');
        
        % bring the images back down to the original size
        LI_1 = LI_1(1:img_h, 1:img_w, :);
        BM_1 = BM_1(1:img_h, 1:img_w, :);
        LI_2 = LI_2(1:img_h, 1:img_w, :);
        BM_2 = BM_2(1:img_h, 1:img_w, :);
        
        % blending the current layer image and the previous image
        img1 = (LI_1.*(1-BM_1) + (img1.*BM_1));
        img2 = (LI_2.*(1-BM_2) + (img2.*BM_2));
        
        bp = 1;
        
    end

    img1 = img1(img_h_range, img_w_range, :);
    img2 = img2(img_h_range, img_w_range, :);
    dm = dm(img_h_range, img_w_range, :);
    
    image_num = num2str(kdx+img_offset, '%03d');
    dm_filename = strcat('depth_maps/dm_', image_num, '.png');
    imwrite(dm, strcat(save_path, dm_filename));

    % now that the image has been created let's run through the intensity values   
    for jdx=1:numel(int_values)
        
        img_int1 = img1*int_values(jdx);
        img_int2 = img2*int_values(jdx);
        
        % save the image file and depth maps
        image_int = num2str(int_values(jdx)*100, '%03d');
        
%         img_filename1 = strcat('images/image_f1_', image_num, '_', image_int, '.png');
%         img_filename2 = strcat('images/image_f2_', image_num, '_', image_int, '.png');
        img_filename1 = strcat('images/image_f1_', image_num, '.png');
        img_filename2 = strcat('images/image_f2_', image_num, '.png');
        
        imwrite(img_int1, strcat(save_path, img_filename1));
        imwrite(img_int2, strcat(save_path, img_filename2));
        
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

