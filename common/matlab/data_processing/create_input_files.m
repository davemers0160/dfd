format long g
format compact
clc
close all
clearvars
%% Info
% use this file to create the text input files for the dnn training
% file contains the following:
%   scenario
%   training data file
%   test data file
%   number of crops
%   train/eval crop sizes
%   input scales
%   number of filters per layer


%% set up global variables
focus_step = 141;
defocus_step = 1;
defocus_start = 128;
defocus_stop = 134;

data_name = 'rw4';
scenario = '6x';
file_suffix = 'left_6x';
cam_side = 'left';

% normal
%version = '_';

% v3
version = '_v3_';

net_version = 'v14a';
stopping_criteria = [60, 400000];
input_root = '/dfd_dnn_rw/inputs';
platform = '';
crop_num = 56;
%crop_size = [32,32, 720,832];
crop_size = [32,32, 624,624];
scale = [1,1];
%con = [256, 64,16,16,64, 256, 128,32,32,128, 256, 256,32,32,256, 128,32,32,128, 64,16,16,64];  %v8d
%con = [256, 128,128,128, 128, 128,128,128, 128, 256,256,256, 256, ...
%            256,256,256, 256, 256,256,256, 256, 256,256,256, 256, ...
%            256,256,256, 256, 128,128,128, 128, 128,128,128]; %v19b
%con = [256, 256,256, 256,256,256, 256, 128,128,128, 128, 128,128,128, 128, 64,64,64]; %v06a
%con = [256, 64,64,64, 64, 64,64,64, 128, 128,128,128, 128, 64,64,64, 64, 64,64,64]; %rw_v17b
con = [256,128,128,128,128,256,256,256,256,512,512,512,256,256,256,128,128,128]; % v14a

%% Get the save folder

save_path = 'D:\IUPUI\DfD\dfd_dnn_rw\';
save_path = uigetdir(save_path, 'Select Input File Save Location');
if(save_path == 0)
    return;
end
commandwindow;

%% Cycle through the voltage steps
for d=defocus_start:defocus_step:defocus_stop
    
    if(d ~= focus_step)
        
        input_file = strcat('dfd_', scenario,'_', cam_side,'_', num2str(focus_step,'%03d_'), num2str(d,'%03d'), version, 'input.txt');
        
        fprintf('Saving: %s\n', input_file);
        file_id = fopen(fullfile(save_path, input_file),'w');
        
        % write the scenario name
        fprintf(file_id, '# Scenario Name\n');
        fprintf(file_id, 'dfd_%s_%s_%03d_%03d%s%02dx%02d_\n', scenario,net_version, focus_step, d, version, crop_size(1), crop_size(2));
        
        % write the stopping criteria
        fprintf(file_id, '# Stopping criteria\n');
        fprintf(file_id, '%d,%d\n', stopping_criteria(1), stopping_criteria(2));
        
        % write the training input: /home/owner/DfD/dfd_dnn_rw/inputs/rw2/d3_8be/dfd_rw2_133_128_train_8b3e.txt
        fprintf(file_id, '# Training/Testing input files\n');
        fprintf(file_id, '%s/%s/%s/dfd_%s_%03d_%03d%strain_%s.txt\n', input_root, data_name, scenario, data_name, focus_step, d, version, file_suffix);
        
        % write the test input: 
        fprintf(file_id, '%s/%s/%s/dfd_%s_%03d_%03d%stest_%s.txt\n', input_root, data_name, scenario, data_name, focus_step, d, version, file_suffix);
        
        % write the number of crops
        fprintf(file_id, '# Number of crops\n');
        fprintf(file_id, '%d\n', crop_num);
        
        % write the training and eval crop sizes
        fprintf(file_id, '# Crop sizes: [training], [evaluation]\n');
        fprintf(file_id, '%d,%d, %d,%d\n', crop_size(1), crop_size(2), crop_size(3), crop_size(4));
        
        % write the scale
        fprintf(file_id, '# Cropping scales/expansion factor\n');
        if(numel(scale)==3)
            fprintf(file_id, '%d,%d, %d\n', scale(1), scale(2), scale(3));
        elseif (numel(scale)==2)
            fprintf(file_id, '%d, %d\n', scale(1), scale(2));
        else
            fprintf(file_id, '%d, %d\n', scale(1), scale(2));
        end
        
        % write the con inputs
        fprintf(file_id, '# %s\n', net_version);
        for idx=1:numel(con)-1
            fprintf(file_id, '%d,', con(idx));
        end
        fprintf(file_id, '%d\n', con(end));
        fclose(file_id);
    end
    
end

fprintf('Complete!\n');




