format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% get the input for the number of blur pixels

b1 = 10;
b2 = 31;



% base image
base_img = cat(1, zeros(100,100), 255*ones(100,100));

box_blur = (1/3)*[1; 1; 1];

figure(plot_num)
% image(base_img(51:250,:))
image(base_img)
colormap(gray(256))
plot_num = plot_num + 1;

%% run each through the blur

b1_img = base_img;
b2_img = base_img;

while(b1_img(ceil(100-b1/2), 100) < 1)

    b1_img = (conv2(b1_img, box_blur, 'same'));
    
end

b1_img = floor(b1_img(51:250, 51:250));

while(b2_img(ceil(100-b2/2), 100) < 1)

    b2_img = (conv2(b2_img, box_blur, 'same'));
    
end

b2_img = floor(b2_img(51:250, 51:250));

figure(plot_num)
hold on
box on
grid on
plot(b1_img(:,100), 'b');
plot(b2_img(:,100), 'g');
plot_num = plot_num + 1;

%%
b1_line = b1_img(:,100) > 0;
b2_line = b2_img(:,100) > 0;

figure(plot_num)
hold on
box on
grid on
plot(b1_line, 'b');
plot(b2_line, 'g');
plot(abs(b1_line - b2_line), 'k')
plot_num = plot_num + 1;

%% image difference

b_diff = abs(b1_img - b2_img);

figure(plot_num)
image(b_diff)
colormap(gray(256))
plot_num = plot_num + 1;

figure(plot_num)
hold on
box on
grid on
plot(b_diff(:,100), 'b');
plot_num = plot_num + 1;






