format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% build the sample data and the blur kernel

% bluring kernel - box filter
% kernel_size = 3;
% kernel = ones(1,kernel_size);
% kernel = kernel/(numel(kernel));

% gaussian kernel size
sigma_start = 0.000;
sigma_step = 0.005;
sigma_stop = 20.0;

%sigma = sigma_start:sigma_step:sigma_stop;

max_blur_radius = 64;
kernel_size = max_blur_radius + 20;

% makesure the kernel size is odd
if(mod(kernel_size,2) == 0)
    kernel_size = kernel_size + 1;
end

commandwindow;

%% run the bluring algorithm

threshold = 1/255;

blur_radius = [0:1:max_blur_radius];

% create a single knife edge line
% 1-D case
data = cat(2, ones(1, 300), 254*ones(1, 600));
% 2-D case
%data = cat(2, zeros(300, 200), 255*ones(300, 400));

fprintf('float[,] kernel = new float[,] {\n');

sigma = sigma_start;
sig_array = zeros(1,numel(blur_radius));

for idx=1:numel(blur_radius)
    
    num = -1;
    while((num < blur_radius(idx)) || (sigma > sigma_stop))
        sigma = sigma + sigma_step;
        
        % 1-D case
        kernel = create_1D_gauss_kernel(kernel_size, sigma);
        blur_data = conv(data, kernel,'same');
        blur_data = double(uint8(blur_data(1:200)));

        % 2-D case
%         kernel = single(create_gauss_kernel(kernel_size, sigma));
%         blur_data = conv2(data, kernel,'same');
%         blur_data = double(floor(blur_data(150:150, 1:400)));
        

        match = (blur_data > (blur_data(1)+threshold)) == (blur_data < (blur_data(end)-threshold));

        num = sum(match);
        
    end
    
    fprintf('{');
    str = '';
    for jdx=floor(kernel_size/2+1):kernel_size
        % 1-D case
        str = strcat(str, num2str(kernel(jdx), '%1.8ff, '));
        % 2-D case
        %str = strcat(str, num2str(kernel(floor(kernel_size/2+1), jdx), '%11.10ff, '));
    end
    str = strcat(str(1:end-1),'},');
    fprintf('%s\n', str);
    
    sig_array(idx) = sigma;
    
end

fprintf('};\n');

%%
fprintf('sigma = [');
str = '';
for idx = 1:numel(sig_array)
    str = strcat(str, num2str(sig_array(idx), '%6.4f, '));
end
str = strcat(str(1:end-1),'];');
fprintf('%s\n', str);

%% 

s2 = sig_array(1:2:end);
s2(1) = 0;

x = [0:2:max_blur_radius];

P = polyfit(x, s2, 2)

slope = x(:)\s2(:)

%% plot the results

% figure(plot_num)
% plot([1:kernel_size], kernel(floor(kernel_size/2),:), 'b')
% hold on
% plot_num = plot_num + 1;

% figure(plot_num)
% plot(data(1:200), 'b')
% hold on
% plot(blur_data, '--r')
% plot_num = plot_num + 1;

figure(plot_num)
plot(sig_array, '.-b')
hold on
plot_num = plot_num + 1;

figure(plot_num)
image(blur_data);
colormap(gray(256));
plot_num = plot_num + 1;

return;

%% This is to create a check of the blur kernel sigma values

max_blur_radius = 64;
kernel_size = 79;

save_path = 'D:\IUPUI\Test_Data\blur_tests';

test_img = cat(2, zeros(300, 100), 255*ones(300, 100+100));

x = 0:1:max_blur_radius;
s3 = 0.0001389*(x.*x) + 0.164*x;

for idx = 1:numel(s3)
    
    kernel = single(create_gauss_kernel(kernel_size, s3(idx)));
    
    blur_img = conv2(test_img, kernel,'same');
    blur_img = blur_img(100:200, 1:200);
    
    imwrite(uint8(blur_img), fullfile(save_path,strcat('blur_image_', num2str(s3(idx), '%07.4f.png'))));
    
end

%% test to display kernel

% =0.0001377901*(A2*A2)+0.163250668*A2
% 0.0001389*x^2 + 0.164*x
max_blur_radius = 64;
kernel_size = 79;

x = 0:1:max_blur_radius;
s3 = 0.0001389*(x.*x) + 0.164*x;

print_half_kernel(s3, kernel_size, 1);

%%
fprintf('sigma = [');
str = '';
for idx = 1:numel(sig_array)
    str = strcat(str, num2str(s3(idx), '%6.4f, '));
end
str = strcat(str(1:end-1),'];');
fprintf('%s\n', str);

