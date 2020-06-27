format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% Data Processing Main Code
% This code will ask the user to select a config file that was generated
% by the GUI to quickly parse all of the images in the raw data

% select a config file - the name of the file will dictate which scenario
% and camera side will be parsed

%load('d:\Common\matlab\lidar\OS1_991827000195.mat');
file_filter = {'*.mat';'*.*'};
data_path = 'D:\Common\matlab\lidar';
[lidar_file, data_path] = uigetfile(file_filter, 'Select LIDAR data file', data_path);
if(data_path == 0)
    return;
end

load(fullfile(data_path,lidar_file));


% get the path to save the results
startpath = 'D:\IUPUI\Test_Data\rw';
save_path = uigetdir(startpath,'Select Save Folder');

if(save_path == 0)
    return;
end

file_filter = {'*.txt','Text Files';'*.*','All Files' };

startpath = 'D:\IUPUI\Test_Data\real_world_raw\configs';
%[cfg_file, cfg_file_path] = uigetfile(file_filter, 'Select Configuration File', startpath);
cfg_file_path = uigetdir(startpath,'Select Config Folder');
if(cfg_file_path == 0)
    return;
end

%% get the listing
v_step = [142,144];
cfg_listing = dir(strcat(cfg_file_path,filesep,'*.txt'));

for idx=1:length(cfg_listing)
    
    process_data(cfg_listing(idx).name, cfg_listing(idx).folder, lidar_struct, v_step, save_path);
    fprintf('------------------------------------------------------\n');
end

fprintf('Complete!\n');

