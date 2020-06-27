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

% http://alumni.media.mit.edu/~wad/color/numbers.html
% color = {[0, 0, 0];[87, 87, 87]/255;[173, 35, 35]/255;[42, 75, 215]/255;...
%          [29, 105, 20]/255;[129, 74, 25]/255;[129, 38, 192]/255;[160, 160, 160]/255;...
%          [129, 197, 122]/255;[157, 175, 255]/255;[41, 208, 208]/255;[255, 146, 51]/255;...
%          [255, 238, 51]/255;[233, 222, 187]/255;[255, 205, 243]/255;[255, 255, 255]/255};
%      
% 6-7-6 RGB color palette https://en.wikipedia.org/wiki/List_of_software_palettes
green = [0, 42, 85, 128, 170, 212, 255];
color = {};
for r=0:5
    for g=0:6
        for b=0:5
            color{end+1} = [51*r, green(g+1), 52*b]/255;
        end
    end
end

commandwindow;

%% create the folders
save_path = 'D:/IUPUI/Test_data/test_blur5/';

warning('off');
mkdir(save_path);
mkdir(save_path, 'images');
mkdir(save_path, 'depth_maps');
warning('on');

%% show the colors
img_w = 50;
img_h = 250;
img = [];

for idx=1:numel(color)
    img = cat(2, img, cat(3, color{idx}(1).*ones(img_h, img_w), color{idx}(2).*ones(img_h, img_w), color{idx}(3).*ones(img_h, img_w))); 
end

figure(plot_num)
imshow(img)
plot_num = plot_num + 1;

%% start to create the images
img_w = 360;
img_h = 360;

% x_min,x_max; y_min,y_max; min_r,max_r
rect = [1,img_w; 1,img_h; 10,30];
circle = [1,img_w; 1,img_h; 10,20];
polygon = [1,img_w; 1,img_h; -50,50];

dm_values = [0, 9:1:232];

fprintf('%s\n\n', save_path);

tic;
for kdx=0:499
    
    [img, dm] = gen_rand_image(img_h, img_w, dm_values, color, {circle, polygon, rect});

%     % get the random background color
%     bg_color = randi([1,numel(color)],1);
%     img = cat(3, color{bg_color}(1).*ones(img_h, img_w), color{bg_color}(2).*ones(img_h, img_w), color{bg_color}(3).*ones(img_h, img_w));
%     dm = zeros(img_h, img_w, 3);
%     
%     D = randi([1, numel(dm_values)],1,50);
%     D = sort(unique(D));
% 
%     for idx=1:numel(D)
% 
%         % get the number of shapes for a given depth map value
% %         N = randi([3,floor((numel(D)-idx)/2) + 6], 1);
%         N = randi([3,floor(exp((4*(numel(D)-idx))/numel(D))) + 6], 1);
% 
%         dm_val = (dm_values(D(idx))/255)*ones(1,3); %, dm_values((idx, idx]/255.0;
%         
%         for jdx=1:N
%             % get the shape type
%             T = randi([1,3], 1);
%             
%             % get the random color for the shape
%             C = randi([1,numel(color)],1);
%             while(C == bg_color)
%                 C = randi([1,numel(color)],1);
%             end
%             
%             switch(T)
%                 case 1
%                     X = randi(circle(1,:), 1);
%                     Y = randi(circle(2,:), 1);
%                     R = randi(circle(3,:), 1);
% 
%                     img = insertShape(img, 'FilledCircle', [X, Y, R], 'Color', color{C}, 'Opacity',1, 'SmoothEdges', false);
%                     dm = insertShape(dm, 'FilledCircle', [X, Y, R], 'Color', dm_val, 'Opacity',1, 'SmoothEdges', false);
% 
%                 case 2
%                     X = randi(polygon(1,:), 1);
%                     Y = randi(polygon(2,:), 1);            
%                     P = randi(polygon(3,:), [1,6]);
%                     P(1:2:end) = P(1:2:end) + X;
%                     P(2:2:end) = P(2:2:end) + Y;
%                     P = cat(2, X, Y, P);
% 
%                     img = insertShape(img, 'FilledPolygon', P, 'Color', color{C}, 'Opacity',1, 'SmoothEdges', false);
%                     dm = insertShape(dm, 'FilledPolygon', P, 'Color', dm_val, 'Opacity',1, 'SmoothEdges', false);
%                     
%                 case 3
%                     X = randi(rect(1,:), 1);
%                     Y = randi(rect(2,:), 1);                      
%                     W = randi(polygon(3,:), 1);
%                     H = randi(polygon(3,:), 1);
%                     
%                     img = insertShape(img, 'FilledRectangle', [X, Y, W, H], 'Color', color{C}, 'Opacity',1, 'SmoothEdges', false);
%                     dm = insertShape(dm, 'FilledRectangle', [X, Y, W, H], 'Color', dm_val, 'Opacity',1, 'SmoothEdges', false);
%             end
%         end
%     end

%     figure(plot_num)
%     imshow(img);
% 
%     figure(plot_num + 1)
%     imshow(dm);

    % save the image file and depth maps
    image_num = num2str(kdx, '%03d');

    img_filename = strcat('images/image_', image_num, '.png');
    imwrite(img, strcat(save_path, img_filename));

    dm_filename = strcat('depth_maps/dm_', image_num, '.png');
    imwrite(dm, strcat(save_path, dm_filename));
    
    fprintf('%s, %s, 0.32, 0.01, 256\n', img_filename, dm_filename);
    
end

toc;


