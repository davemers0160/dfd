format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;


%% select the images

start_path = 'D:\IUPUI\PhD\Results\dfd_dnn_rw\dfd_dnn_tests\';

file_filter = {'*.jpg','JPG Image Files'; '*.png','PNG Image Files'; '*.*','All Files' };

[image_files, image_path] = uigetfile(file_filter, 'Select DFD-Net Log File', start_path, 'multiselect','on');
if(image_path == 0)
    return;
end

save_path = image_path;

commandwindow;

%% start merging images

base_name = 'dfd_rw_sb_k01_k03_merge_';

image_files = sort(image_files);

h_space = 10;

index = 1;

for idx=1:2:numel(image_files)
   
    img1 = imread(fullfile(image_path, image_files{idx}));
    im1_size = size(img1);
    
    % check for an odd number of images to process
    if((idx+1) <= numel(image_files))
        img2 = imread(fullfile(image_path, image_files{idx+1}));
        im2_size = size(img2);
    else
        im2_size = im1_size;
        im2_size(2) = 1;
        img2 = 255*ones(im2_size);
    end
        
    % make sure the images are the same height
    if(im1_size(1) > im2_size(1))
        r_diff = abs(im1_size(1) - im2_size(1));
        img2 = cat(1, img2, 255*ones(r_diff, im2_size(2), im2_size(3)));
    elseif(im1_size(1) < im2_size(1))
        r_diff = abs(im1_size(1) - im2_size(1));
        img1 = cat(1, img1, 255*ones(r_diff, im1_size(2), im1_size(3)));  
    end
    
    pad = 255*ones(im1_size(1), h_space, im1_size(3));
       
    full_img = [img1 pad img2];
       
    save_file = strcat(base_name, 'pt', num2str(index, '%02d'), '.jpg');
    imwrite(full_img, fullfile(save_path, save_file));
    
    % while we're at it let just create the figures for stupid latex
    fprintf('%%----------------------------------------------------------------------\n');
    fprintf('\\begin{sidewaysfigure}[p] \\centering\n');
    fprintf('\\captionsetup[subfloat]{farskip=1pt,captionskip=-14pt}\n');
    fprintf('\\includegraphics[width=8.7in,height=4.7in]{figures/%s} \\\\ [-2mm]\n',save_file);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02da}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02db}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02dc}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02dd}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02de}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02df}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02dg}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\subfloat[\\label{sf:app_dfd_rw_sb_res%02dh}{}]{\\hspace{0.125\\columnwidth}}\n', index);
    fprintf('\\caption{Performance Results for the DfD-Net on the Synthetically Blurred Real World Dataset - Part %d. ', index);
    fprintf('(a) \\& (e) In-focus Image, (b) \\& (f) Out-of-focus Image, (c) \\& (g) Ground Truth Depth Map and (d) \\& (h) DfD-Net Computed Depth Map.}\n');
    fprintf('\\label{fig:app_dfd_rw_sb_res%02d}\n', index);
    fprintf('\\end{sidewaysfigure}\n');
    fprintf('%%----------------------------------------------------------------------\n\n');
     
    index = index + 1;
    
end


