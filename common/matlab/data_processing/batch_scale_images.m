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

startpath = 'D:\IUPUI\Test_Data\rw3';
img_path = uigetdir(startpath,'Select Scenario Folder');

if(img_path == 0)
    return;
end

listing = dir(img_path);
listing = listing(3:end);       % remove .  and .. from the listing

save_dir = img_path;
commandwindow;

%% start the conversion process

img_w = 300;
img_h = 240;

scale_name = '_d3';
fprintf('Path: %s\n',img_path);

for idx=1:length(listing)
    
    if(listing(idx).isdir && (strcmp(listing(idx).name,'left') || strcmp(listing(idx).name,'right')))
        
%         scenario_name = strsplit(listing(idx).folder,filesep);
%         scenario_name = scenario_name{end};
        
        %  create the new folder to put the scaled data: ex left_3
        warning('off','all')
        mkdir(strcat(listing(idx).folder, filesep, listing(idx).name, scale_name));
        warning('on','all')
        
        % get the exposure listing
        exp_listing = dir(strcat(listing(idx).folder,filesep,listing(idx).name,filesep,'exp*'));

        % scan through the exposure folder to scale the images
        for jdx=1:numel(exp_listing)
            exp_name = exp_listing(jdx).name;
            fprintf('Scaling images in %s...\n',exp_name);
            save_dir = strcat(listing(idx).folder, filesep, listing(idx).name, scale_name, filesep, exp_name);
            mkdir(save_dir);
            
            img_listing = dir(strcat(listing(idx).folder, filesep, listing(idx).name, filesep,exp_name, '/*.png'));
            
            for kdx=1:numel(img_listing)
                
                img = imread(fullfile(img_listing(kdx).folder,img_listing(kdx).name));
                img2 = imresize(img,[img_h, img_w], 'lanczos3', 'Antialiasing', true);
                imwrite(img2,fullfile(save_dir,img_listing(kdx).name));
                
            end
            
        end

    end 
    
end

fprintf('Complete!\n');

