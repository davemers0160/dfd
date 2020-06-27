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

load('d:\Common\matlab\lidar\OS1_991827000195.mat');

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

cfg_listing = dir(strcat(cfg_file_path,filesep,'*.txt'));

for idx=1:length(cfg_listing)
    
    cfg_parts = strsplit(cfg_listing(idx).name,{'_','.'});

    scenario_name = cfg_parts{2};
    camera_side = cfg_parts{3};

    folders = strsplit(cfg_file_path,filesep);
    data_path = fullfile(folders{1:end-1},scenario_name);



    [match_params] = parse_input_parameters(fullfile(cfg_file_path, cfg_listing(idx).name));
    process_lidar_data(lidar_struct, match_params, scenario_name, camera_side, data_path, save_path)
    %process_data(cfg_listing(idx).name, cfg_listing(idx).folder, lidar_struct, save_path);
    fprintf('------------------------------------------------------\n');
end

fprintf('Complete!\n');