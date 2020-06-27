format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

%% 

startpath = 'D:\IUPUI\Test_Data\';
%img_path = uigetdir(startpath,'Select Folder Containing Scenarios');
img_path = uigetfolder(startpath, 'Select Folders Containing Scenarios');

if(isempty(img_path))
    return;
end

commandwindow;

%% scan through the directories and get the lidar data

ld_max = 255;
ld_ext = 'left_00000_8bit_x6.png';
ld_search = strcat('/lidar/*',ld_ext);

for idx=1:numel(img_path)
    
    ld_listing = dir(strcat(img_path{idx},'/',ld_search));
    
    if(~isempty(ld_listing))
        [~, ld_name, ~] = fileparts(ld_listing.name);

        ld = imread(fullfile(ld_listing.folder,ld_listing.name));    
        ld_inv = ld_max - ld;

        fprintf('Saving: %s\n',fullfile(ld_listing.folder,strcat(ld_name,'_inv.png')));
        imwrite(ld_inv,fullfile(ld_listing.folder,strcat(ld_name,'_inv.png')));
    end
end