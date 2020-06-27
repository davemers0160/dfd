format long g
format compact
clc
%close all
%clearvars

%% got through and create the image inputs
save_path = 'D:\IUPUI\DfD\dfd_dnn_rw\inputs\test\';

start_step = 140;
end_step = 129;

% this is a list of the scenarios to comment out for training and testing
%training
%comment = {'Auditorium2','GR2','House1','IWM3','WS2'};
%testing
comment = {'Auditorium1','GR2','Garage1','Garage2','Garage3','Grass1','Grass2','IWM1','IWM2','IWM4','IWM5','IWM6','IWM7','Lab1','Lab2','Library1','Library2','Library3','SL012','Shed1','Shed2','WS1','WS3','WS4'};

%% get the data listing for the scenarios
startpath = 'D:\IUPUI\Test_Data\rw';

data_path = uigetdir(startpath,'Select Data Folder');

if(data_path == 0)
    return;
end

listing = dir(data_path);
listing = listing(3:end);

%%
for idx=start_step:-1:end_step+1
    for jdx=end_step:1:idx-1
        
        file_name = strcat('dfd_rw_',num2str(idx),'_',num2str(jdx),'_test_input_v2.txt');
        file_id = fopen(fullfile(save_path, file_name),'w');
        
        fprintf(file_id, '# Full data listing for the real world dataset\n# Data Directory: \n');
        fprintf(file_id, 'D:/IUPUI/Test_Data/rw/, /home/owner/DfD/data/rw/, /home/owner/DfD/data/rw/\n');
        fprintf(file_id, '# Data: \n');        
        for kdx=1:length(listing)
            
            if(listing(kdx).isdir)
                scenario_name = listing(kdx).name;
        
                if(any(strcmp(scenario_name,comment)))
                    c='#';
                else
                    c=''; 
                end
                % get the left exposure level listing
                left_exp_listing = dir(strcat(fullfile(listing(kdx).folder,listing(kdx).name,'left'),filesep,'exp*'));

                for mdx=1:length(left_exp_listing)
                    exp_name = left_exp_listing(mdx).name;
                    exp_split = strsplit(exp_name,'_');
                    exp_split = str2double(exp_split{2});
                    focus_file = strcat(scenario_name,'/left/',exp_name,'/image_',num2str(idx),'_',num2str(exp_split),'.00.png');
                    defocus_file = strcat(scenario_name,'/left/',exp_name,'/image_',num2str(jdx),'_',num2str(exp_split),'.00.png');
                    fprintf(file_id, '%s%s, %s, %s\n', c, focus_file, defocus_file, strcat(scenario_name,'/lidar/lidar_rng_left_00000_8bit.png'));
                end

                fprintf(file_id, '#------------------------------------------------------\n');
                
                % get the right exposure level listing
                right_exp_listing = dir(strcat(fullfile(listing(kdx).folder,listing(kdx).name,'right'),filesep,'exp*'));

                for mdx=1:length(right_exp_listing)
                    exp_name = right_exp_listing(mdx).name;
                    exp_split = strsplit(exp_name,'_');
                    exp_split = str2double(exp_split{2});
                    focus_file = strcat(scenario_name,'/right/',exp_name,'/image_',num2str(idx),'_',num2str(exp_split),'.00.png');
                    defocus_file = strcat(scenario_name,'/right/',exp_name,'/image_',num2str(jdx),'_',num2str(exp_split),'.00.png');
                    fprintf(file_id, '%s%s, %s, %s\n', c, focus_file, defocus_file, strcat(scenario_name,'/lidar/lidar_rng_right_00000_8bit.png'));
                end
                
                fprintf(file_id, '#------------------------------------------------------\n');

            end
          
        end
        fclose(file_id);

    end    
end
fprintf('Complete!\n');


%% create the listing of the pairs


for idx=start_step:-1:end_step+1
    for jdx=end_step:1:idx-1        
        fprintf('{''dfd_rw_%d_%d_train_input_v2.txt'',''dfd_rw_%d_%d_test_input_v2.txt''};\n',idx, jdx, idx, jdx);
    end    
end