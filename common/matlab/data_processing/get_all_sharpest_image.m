%% this function looks at a scenario and picks the most infocus image

format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

%% Select all of the folders

startpath = 'D:\IUPUI\Test_Data\'; 
[folder_names] = uigetfolder(startpath, 'Select Scenario Folders to compute Sharpest Image');

camera_side = {'left','right'};

commandwindow;

%% cycle through the selected directories

for idx=1:numel(folder_names)
    
    folders = strsplit(folder_names{idx}, filesep);
    scenario_name = folders{end};

    get_sharpest_image(folder_names{idx}, folders, scenario_name, camera_side);
    
end

fprintf('Complete!\n');
