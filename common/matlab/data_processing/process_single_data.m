format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

%% load lidar struct
%load('d:\Common\matlab\lidar\OS1_991827000195.mat');
% file_filter = {'*.mat','MAT Files';'*.*','All Files'};
% lidar_data_path = 'D:\Common\matlab\lidar';
% [lidar_file, lidar_data_path] = uigetfile(file_filter, 'Select LIDAR data file', lidar_data_path);
% if(lidar_data_path == 0)
%     return;
% end

% load(fullfile(lidar_data_path,lidar_file));

%% Load config files
% get the path to save the results
startpath = 'D:\IUPUI\Test_Data\';
save_path = uigetdir(startpath,'Select Save Folder');

if(save_path == 0)
    return;
end

file_filter = {'*.txt','Text Files';'*.*','All Files' };

[cfg_file, cfg_file_path] = uigetfile(file_filter, 'Select Configuration File', startpath, 'MultiSelect', 'on');
if(cfg_file_path == 0)
    return;
end
v_step = [122,148];
commandwindow;

if(iscell(cfg_file))
    for idx=1:numel(cfg_file)
        process_data(cfg_file{idx}, cfg_file_path(1:end-1), v_step, save_path);
        fprintf('---------------------------------------------------------\n');
    end
else
    process_data(cfg_file, cfg_file_path(1:end-1), v_step, save_path);
end
fprintf('Complete!\n');

