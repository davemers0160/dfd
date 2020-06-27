format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

%% Process a single file
%load('d:\Common\matlab\lidar\OS1_991827000195.mat');
file_filter = {'*.mat','MAT Files';'*.*','All Files'};
lidar_data_path = 'D:\Common\matlab\lidar';
[lidar_file, lidar_data_path] = uigetfile(file_filter, 'Select LIDAR data file', lidar_data_path);
if(lidar_data_path == 0)
    return;
end

load(fullfile(lidar_data_path,lidar_file));

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
commandwindow;

%% parse the contents of the file

if(iscell(cfg_file))
    
    for idx=1:numel(cfg_file)
        cfg_parts = strsplit(cfg_file{idx},{'_','.'});

        scenario_name = cfg_parts{2};
        camera_side = cfg_parts{3};

        folders = strsplit(cfg_file_path,filesep);
        data_path = fullfile(folders{1:end-2},scenario_name);

        [match_params] = parse_input_parameters(fullfile(cfg_file_path, cfg_file{idx}));
        process_lidar_data(lidar_struct, match_params, scenario_name, camera_side, data_path, save_path)
        fprintf('---------------------------------------------------------\n');
    end

else
    cfg_parts = strsplit(cfg_file,{'_','.'});

    scenario_name = cfg_parts{2};
    camera_side = cfg_parts{3};

    folders = strsplit(cfg_file_path,filesep);
    data_path = fullfile(folders{1:end-2},scenario_name);

    [match_params] = parse_input_parameters(fullfile(cfg_file_path, cfg_file));
    process_lidar_data(lidar_struct, match_params, scenario_name, camera_side, data_path, save_path)

end
fprintf('Complete!\n');