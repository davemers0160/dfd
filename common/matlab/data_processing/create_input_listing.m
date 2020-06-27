format long g
format compact
clc
close all
clearvars
%% info
% use this file to create the actual data input files for testing and
% training

%% set up global variables
focus_step = 141;
defocus_step = 1;
defocus_start = 128;    %124
defocus_stop = 134;     %147
defocus_start2 = [];    %125;

cam_side_left = 'left';
cam_side_right = 'none';
lidar_bit = '8bit_x6';
data_name = 'rw4';
file_suffix = 'left_6x';
% normal
%train_version = 'train_';
%test_version = 'test_';
% v3
train_version = 'v3_train_';
test_version = 'v3_test_';


%% Get the save folder

save_path = 'D:\IUPUI\DfD\dfd_dnn_rw\inputs\';
save_path = uigetdir(save_path, 'Select Input File Save Location');
if(save_path == 0)
    return;
end

%% Select the folders containing the training data
start_path = 'D:\IUPUI\Test_Data\';
training_data = uigetfolder(start_path, 'Select Training Data Folders');

if(isempty(training_data))
    return;
end
%% Select the folders containing the test data
start_path = 'D:\IUPUI\Test_Data\';
test_data = uigetfolder(start_path, 'Select Test Data Folders');

if(isempty(test_data))
    return;
end
%% Cycle through the voltage steps
commandwindow;
for d=defocus_start:defocus_step:defocus_stop
    
    if(d ~= focus_step)
        % set up the file save names format
        training_file = strcat('dfd_', data_name, '_', num2str(focus_step,'%03d_'), num2str(d,'%03d_'), train_version, file_suffix, '.txt');
        test_file = strcat('dfd_', data_name, '_', num2str(focus_step,'%03d_'), num2str(d,'%03d_'), test_version, file_suffix, '.txt');

        fprintf('Saving: %s\n',training_file);
        file_id = fopen(fullfile(save_path, training_file),'w');

        fprintf(file_id, '# Full data listing for the real world dataset\n# Data Directory: \n');
        fprintf(file_id, 'Test_Data/%s/, DfD/data/%s/, Projects/rw_data/\n', data_name, data_name);
        fprintf(file_id, '# Data: \n');   
        fprintf(file_id, '#------------------------------------------------------\n');

        % cycle through each folder and save the results
        for idx=1:numel(training_data)

            listing = dir(training_data{idx});
            listing = listing(3:end);

            for jdx=1:numel(listing)

                if(listing(jdx).isdir && (strcmp(listing(jdx).name,cam_side_left) || strcmp(listing(jdx).name,cam_side_right)))
                    % get the scenario name
                    if(strcmp(listing(jdx).name,cam_side_left))
                        cam_side = 'left';
                    else
                        cam_side = 'right';
                    end
                    
                    scenario_name = strsplit(listing(jdx).folder,filesep);
                    scenario_name = scenario_name{end};

                    % get the exposure folder names
                    exp_listing = dir(strcat(training_data{idx},filesep,listing(jdx).name,filesep,'exp*'));

                    % cycle through each exposure folder name and generate an
                    % image tuple
                    for kdx=1:length(exp_listing)
                        exp_name = exp_listing(kdx).name;                
                        exp_parts = strsplit(exp_name, {'_'});

                        focus_file = strcat(scenario_name,'/',listing(jdx).name,'/',exp_name,'/image_',num2str(focus_step),'_',exp_parts{2},'.00.png');
                        defocus_file = strcat(scenario_name,'/',listing(jdx).name,'/',exp_name,'/image_',num2str(d),'_',exp_parts{2},'.00.png');
                        if(~isempty(defocus_start2))
                            defocus_file2 = strcat(scenario_name,'/',listing(jdx).name,'/',exp_name,'/image_',num2str(defocus_start2),'_',exp_parts{2},'.00.png');
                            fprintf(file_id, '%s, %s, %s, %s\n', focus_file, defocus_file, defocus_file2, strcat(scenario_name,'/lidar/lidar_rng_',cam_side,'_00000_',lidar_bit,'.png'));
                        else
                            fprintf(file_id, '%s, %s, %s\n', focus_file, defocus_file, strcat(scenario_name,'/lidar/lidar_rng_',cam_side,'_00000_',lidar_bit,'.png'));                            
                        end
                    end

                    fprintf(file_id, '#------------------------------------------------------\n');

                end
            end

        end
        fclose(file_id);

        fprintf('Saving: %s\n',test_file);
        file_id = fopen(fullfile(save_path, test_file),'w');

        fprintf(file_id, '# Full data listing for the real world dataset\n# Data Directory: \n');
        fprintf(file_id, 'Test_Data/%s/, DfD/data/%s/, Projects/rw_data/\n', data_name, data_name);
        fprintf(file_id, '# Data: \n');   
        fprintf(file_id, '#------------------------------------------------------\n');

        % cycle through each folder and save the results
        for idx=1:numel(test_data)

            listing = dir(test_data{idx});
            listing = listing(3:end);

            for jdx=1:numel(listing)

                if(listing(jdx).isdir && (strcmp(listing(jdx).name, cam_side_left) || strcmp(listing(jdx).name, cam_side_right)))
                    % get the scenario name
                    scenario_name = strsplit(listing(jdx).folder,filesep);
                    scenario_name = scenario_name{end};

                    % get the exposure folder names
                    exp_listing = dir(strcat(test_data{idx},filesep,listing(jdx).name,filesep,'exp*'));

                    % cycle through each exposure folder name and generate an
                    % image tuple
                    for kdx=1:length(exp_listing)
                        exp_name = exp_listing(kdx).name;                
                        exp_parts = strsplit(exp_name, {'_'});

                        focus_file = strcat(scenario_name,'/',listing(jdx).name,'/',exp_name,'/image_',num2str(focus_step),'_',exp_parts{2},'.00.png');
                        defocus_file = strcat(scenario_name,'/',listing(jdx).name,'/',exp_name,'/image_',num2str(d),'_',exp_parts{2},'.00.png');
                        if(~isempty(defocus_start2))
                            defocus_file2 = strcat(scenario_name,'/',listing(jdx).name,'/',exp_name,'/image_',num2str(defocus_start2),'_',exp_parts{2},'.00.png');
                            fprintf(file_id, '%s, %s, %s, %s\n', focus_file, defocus_file, defocus_file2, strcat(scenario_name,'/lidar/lidar_rng_',cam_side,'_00000_',lidar_bit,'.png'));
                        else
                            fprintf(file_id, '%s, %s, %s\n', focus_file, defocus_file, strcat(scenario_name,'/lidar/lidar_rng_',cam_side,'_00000_',lidar_bit,'.png'));
                        end
                    end

                    fprintf(file_id, '#------------------------------------------------------\n');

                end
            end

        end
        fclose(file_id);

    end
    defocus_start2 = defocus_start2 + 2;
end
fprintf('Complete!\n');

