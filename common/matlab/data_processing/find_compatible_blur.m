format long g
format compact
%clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% get the two images

file_filter = {'*.png','PNG Image Files';'*.*','All Files' };

if(ispc)
    startpath = 'D:\IUPUI\Test_Data\rw_raw4\k64\left';
else
    startpath = '/home/owner/DfD';
end

[infocus_img_file, img_path1] = uigetfile(file_filter, 'Select In-Focus Image File', startpath);
if(img_path1 == 0)
    return;
end

[defocus_img_file, img_path2] = uigetfile(file_filter, 'Select Out-of-Focus Image File', img_path1);
if(img_path2 == 0)
    return;
end

commandwindow;

fprintf('%s\n', fullfile(img_path2,defocus_img_file));

%% process the images
kernel_size = 19;
sigma_start = 0.20;
sigma_stop = 3.5;
sigma_step = 0.01;
sigma = sigma_start:sigma_step:sigma_stop;

% set the crop width, use an even number
% crop_w = 100;   
%crop_h = 140:900; %large
crop_h = 19:728;    %small


% read in the images
infocus_img = double(imread(fullfile(img_path1,infocus_img_file)));
defocus_img = double(imread(fullfile(img_path2,defocus_img_file)));

%infocus_img = infocus_img(19:end,:,:);
%defocus_img = defocus_img(19:end,:,:);

infocus_img = infocus_img(crop_h,:,:);
defocus_img = defocus_img(crop_h,:,:);

infocus_img_g = double(rgb2gray(uint8(infocus_img)));
defocus_img_g = double(rgb2gray(uint8(defocus_img)));

figure(plot_num); 
image(infocus_img_g);
colormap(gray(256));
title('In-focus Image');
plot_num = plot_num + 1;

figure(plot_num);
image(defocus_img_g);
colormap(gray(256));
title('Out-of-Focus Image');
plot_num = plot_num + 1;

% get the center section of the in-focus image and out-of-focus image
[img_h, img_w] = size(infocus_img_g);

% [x, y]
%ring_center = [306, 178];    % 306,196 without the 18 removed from the top
ring_center = [360, 313];
%ring_center = [625,200];     % 570,340 for the larger image

r_max = [sqrt((ring_center(1)-1)^2 + (ring_center(2)-1)^2), sqrt((ring_center(1)-1)^2 + (ring_center(2)-img_h)^2), ...
         sqrt((ring_center(1)-img_w)^2 + (ring_center(2)-1)^2), sqrt((ring_center(1)-img_w)^2 + (ring_center(2)-img_h)^2)];

r_step = 20;
r_outer = 20:r_step:ceil(max(r_max)+2*r_step);
r_inner = 0;

% build the mask parameters based on the ring parameters
mask_params.center = ring_center;
mask_params.range = r_outer;
mask_parmas.step = r_step;


%% blur the defocused image some more
% nk = createGaussKernel(kernel_size, 3.0);
% full_blur_img = imfilter(infocus_img_g, nk, 'symmetric', 'conv');
% 
% 
% figure(plot_num);
% image(full_blur_img);
% colormap(gray(256));
% title('Blurred In-Focus Image');
% plot_num = plot_num + 1;

% crop_in = infocus_img(floor(img_h/2-crop_h/2):floor(img_h/2+crop_h/2)-1, floor(img_w/2-crop_w/2):floor(img_w/2+crop_w/2)-1);
% crop_de = defocus_img(floor(img_h/2-crop_h/2):floor(img_h/2+crop_h/2)-1, floor(img_w/2-crop_w/2):floor(img_w/2+crop_w/2)-1);

%% cycle through the possible options for blur filters

mask = create_ring_mask([img_h,img_w], ring_center, r_outer(1), 0);
mask(ring_center(2),ring_center(1)) = 1;    % do this becasue of the math for removing the center radius
kernel = {};
best_mse_index = 91;
best_min_mse = 0;

mask = ones(img_h, img_w);

[full_blur_img, kernel{1}, best_mse_index, best_min_mse] = get_best_blur_kernel(infocus_img_g, defocus_img_g, mask, kernel_size, sigma);
%[full_blur_img, kernel{1}, best_mse_index, best_min_mse] = get_best_blur_kernel(blur_img, defocus_img_g, mask, kernel_size, sigma);
%[full_blur_img, kernel{1}, best_mse_index, best_min_mse] = get_best_blur_kernel(infocus_img_g, blur_img, mask, kernel_size, sigma);

fprintf('best_min_mse[%03d]: %2.4f;    sigma: %2.4f\n', best_mse_index, best_min_mse, sigma(best_mse_index));

figure(plot_num);
image(full_blur_img);
colormap(gray(256));
title('Blurred In-Focus Image');
plot_num = plot_num + 1;


%% blur the entire image with the best matching kernel

% kernel{1} = createGaussKernel(kernel_size, 1.0);
% full_blur_img = imfilter(infocus_img, kernel{1}, 'symmetric', 'conv');

%img_diff = (blur_img - full_blur_img);

m = ceil(kernel_size/2);

img_diff = (defocus_img_g - full_blur_img);
plot_img_diff = abs(img_diff(m:end-m,m:end-m));

img_mse_baseline = mean(img_diff(:).*img_diff(:));

figure(plot_num);
set(gcf,'position',([50,50,900,700]),'color','w', 'Name','Original Image Difference');

%surf(abs(img_diff(m:end-m,m:end-m)));
surf(plot_img_diff);

%surf(abs(img_diff(m:end-m,m:end-m) + 20*mask(m:end-m,m:end-m)));
shading interp
colormap(jet(512));

% x-Axis
xlim([1 size(plot_img_diff,2)]);
xticklabels([]);
xlabel('X', 'fontweight','bold','FontSize', 13);

% y-axis
ylim([1 size(plot_img_diff,1)]);
yticklabels([]);
ylabel('Y', 'fontweight','bold','FontSize', 13);

% z-axis
zlabel('Abs Difference', 'fontweight','bold','FontSize', 13);

set(gca, 'fontweight','bold');
grid on;
box on;

ax = gca;
ax.YDir = 'reverse';
ax.TickLength=[0,0];

%title(strcat('\middefocus\_img - blur\_img\mid; sigma:',32,num2str(sigma(best_mse_index),'%2.3f')),'Interpreter','tex');
view(-30,70);

ax.Position = [0.07 0.03 0.9 0.94];

plot_num = plot_num + 1;

%% plot the same thing just as a contour plot
w=11;
n = 15;
k = ones(w)*(1/(w*w));
tmp_diff = imfilter(plot_img_diff,k,'replicate');

tmp_diff = ceil(tmp_diff*n)/n;

figure(plot_num);
set(gcf,'position',([50,50,900,700]),'color','w', 'Name','Original Image Difference');

% t2 = tmp_diff;
% t2(tmp_diff<5.5) = 200;

mt=min(tmp_diff(:));
%contour(tmp_diff,[0,mt+1]);
surf(tmp_diff);
shading interp;
colormap(jet(n));

zlim([0, 30]);

ax = gca;
ax.YDir = 'reverse';

plot_num = plot_num + 1;


%%

% figure(plot_num);
% set(gcf,'position',([50,50,900,700]),'color','w')
% mask = create_ring_mask([img_h,img_w], [306,196], 150, 0);
% 
% imagesc(img_diff(m:end-m,m:end-m)+20*mask(m:end-m,m:end-m));
% 
% %shading interp
% colormap(jet(512));
% 
% ax = gca;
% %ax.YDir = 'reverse';
% 
% title(defocus_img_file(1:9),'interpreter','none')
% 
% plot_num = plot_num + 1;

%%

if(false)
    img_diff2 = defocus_img_g - infocus_img_g;

    figure(plot_num);
    set(gcf,'position',([50,50,900,700]),'color','w')

    surf(img_diff2(m:end-m,m:end-m)+20*mask(m:end-m,m:end-m));
    shading interp
    colormap(jet(512));

    ax = gca;
    ax.YDir = 'reverse';

    title(defocus_img_file(1:9),'interpreter','none')
    view(-30,70);

    plot_num = plot_num + 1;
end
%% run through the different radius' and find the best blur for each ring

mask = create_ring_mask([img_h,img_w], ring_center, r_outer(1), 0);
mask(ring_center(2),ring_center(1)) = 1;    % do this becasue of the math for removing the center radius
corr_img = [];
corr_img = full_blur_img.*mask;

%tmp = zeros(img_h,img_w);

mse_index(1) = best_mse_index;
min_mse(1) = best_min_mse;
%figure(plot_num);
for idx=2:numel(r_outer)-1
    
    %create the mask
    mask = create_ring_mask([img_h,img_w], ring_center, r_outer(idx), r_outer(idx-1));

    %[b, kernel{idx}, mse_index(idx), min_mse(idx)] = get_best_blur_kernel(infocus_img_g, defocus_img_g, mask, kernel_size, sigma);    
    %[b, kernel{idx}, mse_index(idx), min_mse(idx)] = get_best_blur_kernel(infocus_img_g, blur_img, mask, kernel_size, sigma);    
    [b, kernel{idx}, mse_index(idx), min_mse(idx)] = get_best_blur_kernel(defocus_img_g, full_blur_img, mask, kernel_size, sigma);    

    fprintf('%04d: min_mse[%03d]: %2.4f;    sigma: %2.4f\n', idx, mse_index(idx), min_mse(idx), sigma(mse_index(idx)));

%     tmp_blur_img = imfilter(infocus_img_g, kernel{idx}, 'symmetric', 'conv');
%     %tmp_blur_img = imfilter(blur_img, kernel{idx}, 'symmetric', 'conv');
% 
%     corr_img = corr_img + (tmp_blur_img.*mask);
    %     tmp = tmp+(mask*idx);
    %surf(img_diff(m:end-m,m:end-m)-1*tmp(m:end-m,m:end-m));
%     surf(corr_img);
%     shading interp
%     ax = gca;
%     ax.YDir = 'reverse';
%     view(-30,70);
%     drawnow;
%      pause(0.1);
end

%plot_num = plot_num + 1;

fprintf('Ring Blur Search Complete!\n');

%% build the new image...

mask = create_ring_mask([img_h,img_w], ring_center, r_outer(1), 0);
mask(ring_center(2),ring_center(1)) = 1;    % do this becasue of the math for removing the center radius

%nk = createGaussKernel(kernel_size, sigma(1));
tmp_blur_img = imfilter(defocus_img_g, kernel{1}, 'symmetric', 'conv');

corr_img_g = [];
corr_img_g = tmp_blur_img.*mask;

for idx=2:numel(r_outer)-1
    
%     if(abs(mse_index(1) ~= mse_index(idx)))
%         nk = createGaussKernel(kernel_size, sigma(abs(mse_index(1)-mse_index(idx))));
        tmp_blur_img = imfilter(defocus_img_g, kernel{idx}, 'symmetric', 'conv');
%     else
%         nk = createGaussKernel(kernel_size, sigma(1));
%         tmp_blur_img = imfilter(defocus_img, nk, 'symmetric', 'conv');
%     end
    
    mask = create_ring_mask([img_h,img_w], ring_center, r_outer(idx), r_outer(idx-1));
    corr_img_g = corr_img_g + (tmp_blur_img.*mask);

end


% figure(plot_num);
% image(defocus_img);
% colormap(gray(256));
% title('defocus-img');
% plot_num = plot_num + 1;

figure(plot_num);
image(corr_img_g);
colormap(gray(256));
title('corr-img');
plot_num = plot_num + 1;


%% run a check by redoing the difference between the infocus image and the new corrected defocused image


% mask = create_ring_mask([img_h,img_w], ring_center, r_outer(1), 0);
% mask(ring_center(2),ring_center(1)) = 1;    % do this becasue of the math for removing the center radius
%kernel = {};

%[b, kernel{1}, best_mse_index] = get_best_blur_kernel(infocus_img, corr_img, mask, kernel_size, sigma);
%fprintf('best_min_mse[%d]: sigma: %2.4f\n', best_mse_index, sigma(best_mse_index));

full_blur_img = imfilter(infocus_img_g, kernel{1}, 'symmetric', 'conv');
%full_blur_img = imfilter(blur_img, kernel{1}, 'symmetric', 'conv');

m = ceil(kernel_size/2);

img_diff = (corr_img_g - full_blur_img);
plot_img_diff = abs(img_diff(m:end-m,m:end-m));

img_mse = img_diff.*img_diff;
img_mse = mean(img_mse(:));

fprintf('Baseline Blur Image MSE: %2.4f\n',img_mse_baseline);
fprintf('Corrected Blur Image MSE: %2.4f\n',img_mse);


figure(plot_num);
set(gcf,'position',([50,50,900,700]),'color','w','Name','Corrected Image Difference')

surf(plot_img_diff);
%surf(abs(img_diff));
shading interp;
colormap(jet(512));

% x-Axis
xlim([1 size(plot_img_diff,2)]);
xticklabels([]);
xlabel('X', 'fontweight','bold','FontSize', 13);

% y-axis
ylim([1 size(plot_img_diff,1)]);
yticklabels([]);
ylabel('Y', 'fontweight','bold','FontSize', 13);

% z-axis
zlabel('Abs Difference', 'fontweight','bold','FontSize', 13);

set(gca, 'fontweight','bold');
box on;
grid on;

ax = gca;
ax.YDir = 'reverse';
ax.TickLength=[0,0];

%title(strcat('corr_img - full_blur_img',32,defocus_img_file(1:9)),'interpreter','none');
%title(strcat('\midcorrected\_img - blur\_img\mid; sigma:',32,num2str(sigma(best_mse_index),'%2.3f')),'Interpreter','tex');
view(-30,70);

ax.Position = [0.07 0.03 0.9 0.94];

plot_num = plot_num + 1;

%%
% nk = createGaussKernel(kernel_size, sigma(1));
% new_infocus_img = imfilter(infocus_img, nk, 'symmetric', 'conv');
% new_defocus_img = imfilter(defocus_img, nk, 'symmetric', 'conv');
% 
% figure; image(new_infocus_img);colormap(gray(256));title('new-infocus-img');
% figure; image(new_defocus_img);colormap(gray(256));title('new-defocus-img');
% 
% plot_num = plot_num + 2;
% 
% 
% img_diff = new_infocus_img-new_defocus_img;

% img_diff2 = corr_img - infocus_img_g;
% 
% 
% 
% figure(plot_num);
% set(gcf,'position',([50,50,900,700]),'color','w')
% 
% surf(abs(img_diff2(m:end-m,m:end-m)));
% shading interp
% colormap(jet(512));
% 
% ax = gca;
% ax.YDir = 'reverse';
% 
% title(strcat('img diff',32,defocus_img_file(1:9)),'interpreter','none')
% view(-30,70);
% 
% plot_num = plot_num + 1;

%% save the important data to a mat file
save_vars = false;

if(save_vars)
    mat_save_path = 'D:\IUPUI\Test_Data\rw_raw4\k64';
    mat_save_file = 'k64_128_blur_correction_params.mat';
    save(fullfile(mat_save_path, mat_save_file), 'kernel', 'mask_params');
end

%% try the correction on the color version of the image

[corr_img] = apply_blur_correction(defocus_img, kernel, mask_params);

figure(plot_num);
image(uint8(corr_img));
%colormap(gray(256));
title('corr-img');
plot_num = plot_num + 1;

%%
% figure; image(corr_img);colormap(gray(256));title('corr-img');
% %figure; image(new_blur_img);colormap(gray(256));title('new-blur-img');
% figure; image(defocus_img);colormap(gray(256));title('defocus-img');
% 
% figure; image(infocus_img);colormap(gray(256));title('infocus-img');
% figure; surf(infocus_img-corr_img);colormap(gray(256));shading interp;title('diff-img');



