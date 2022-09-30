format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[scriptpath,  filename, ext] = fileparts(full_path);
plot_count = 1;
line_width = 1.0;

commandwindow;

%%
sigma = [0.1725,0.3450,0.5175,0.6900,0.8625,1.0350,1.2075,1.3800,1.5525,1.7250,1.8975,2.0700,2.2425,2.4150,2.5875,2.7600,2.9325,3.1050,3.2775,3.4500,3.6225,3.7950,3.9675,4.1400,4.3125,4.4850,4.6575,4.8300,5.0025];

num_images = 14;
kernel_size = 69;

% create the knife edge image
img = cat(2, 255*ones(300,150), zeros(300,150));


%% cycle through the sigmas and save the images

for idx=1:num_images
    
    k = create_gauss_kernel(kernel_size, sigma(idx));
    
    img_f = conv2(img, k, 'same');
    
    img_f = uint8(img_f(51:250, 51:250));
    
    image(img_f);
    colormap(gray(256));
    
    imwrite(img_f, strcat('D:/data/dfd/blur_example/blur_example_i', num2str(idx-1, '%02d'),'.png'));
    
    pause(1);
    
end

