%% Process Unity DfD Images
format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% get the input file

file_filter = {'*.txt','Text Files';'*.*','All Files' };
startpath = 'D:\Projects\';

[input_file, input_path] = uigetfile(file_filter, 'Select Data Input File', startpath);
if(input_path == 0)
    return;
end

% get the save directory
save_path = uigetdir(input_path, 'Select Save Folder');

if(save_path == 0)
    return;
end

warning('off');
mkdir(save_path, 'images');
mkdir(save_path, 'depthmaps');
warning('on');

commandwindow;

%% read in the file

params = parse_input_parameters(fullfile(input_path, input_file));

% get the directory for the data
data_directory = params{1}{1};
params(1) = [];

%% run through the images and crop them and then convert the dat file into a png

x_offset = 20;
y_offset = 20;

red = zeros(2,1);
green = zeros(2,1);
blue = zeros(2,1);
gry = zeros(2,1);

str_line = cell(length(params),1);

for idx=1:length(params)    

    % this is expected to be a 3-channel color image
    img1 = imread(fullfile(data_directory, params{idx}{1}));
    img2 = imread(fullfile(data_directory, params{idx}{2}));
    dm = read_binary_depthmap_file(fullfile(data_directory, params{idx}{3}), 'uint32', 'uint16');
    [~,dm_filename,~] = fileparts(params{idx}{3});
    
    str_line{idx} = strcat(params{idx}{1}, ',', 32, params{idx}{2}, ',', 32, 'depthmaps/', dm_filename,'.png');
    
    [img_h, img_w, ~] = size(img1);
       
    % crop the images
    img1 = img1(y_offset:(img_h-y_offset), x_offset:(img_w-x_offset), 1:3);
    img2 = img2(y_offset:(img_h-y_offset), x_offset:(img_w-x_offset), 1:3);
    dm = dm(y_offset:(img_h-y_offset), x_offset:(img_w-x_offset));
    
        % get the combined color stats
    r = mean(mean(img1(:,:,1)));
    g = mean(mean(img1(:,:,2)));
    b = mean(mean(img1(:,:,3)));
    gr = mean(mean(rgb2gray(img1)));
    
    red(1) = red(1) + r;
    green(1) = green(1) + g; 
    blue(1) = blue(1) + b;
    gry(1) = gry(1) + gr;
    
    r = mean(mean(img2(:,:,1)));
    g = mean(mean(img2(:,:,2)));
    b = mean(mean(img2(:,:,3)));
    gr = mean(mean(rgb2gray(img2)));
    
    red(2) = red(2) + r;
    green(2) = green(2) + g; 
    blue(2) = blue(2) + b;
    gry(2) = gry(2) + gr;    

    image(dm);
    colormap(jet(20));
    
%     imwrite(img1, fullfile(save_path, params{idx}{1}));
%     imwrite(img2, fullfile(save_path, params{idx}{2}));
%     imwrite(uint8(dm), fullfile(save_path, 'depthmaps', strcat(dm_filename,'.png')));
%         
    input(' ')
end

red = red/length(params);
green = green/length(params);
blue = blue/length(params);
gry = gry/length(params);

%% write the input file

save_name = strcat(save_path, '/input_file.txt');
file_id = fopen(save_name, 'w');
%fprintf(file_id, '# %s\n\n', color_palette);
%fprintf(file_id, '%s\n\n', data_directory);

fprintf(file_id, '# color stats: r1, g1, b1, r2, g2, b2, gray1, gray2\n');
fprintf(file_id, '# color stats: %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n\n', red(1), green(1), blue(1), red(2), green(2), blue(2), gry(1), gry(2));

fprintf('# color stats: r1, g1, b1, r2, g2, b2, gray1, gray2\n');
fprintf('# color stats: %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n\n', red(1), green(1), blue(1), red(2), green(2), blue(2), gry(1), gry(2));

save_path = strrep(save_path, '\', '/');
if(save_path(end) ~= '/')
    save_path = strcat(save_path, '/');
end

fprintf(file_id, '%s\n\n', save_path);
fprintf('%s\n\n', save_path);

for idx=1:length(params) 
    
    fprintf(file_id, '%s\n', str_line{idx});
    fprintf('%s\n', str_line{idx});
end

fprintf(file_id, '\n');
fprintf('\n');

fclose(file_id);
