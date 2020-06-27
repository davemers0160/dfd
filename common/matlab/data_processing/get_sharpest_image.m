% %% this function looks at a scenario and picks the most infocus image
% 
% format long g
% format compact
% clc
% close all
% clearvars
% 
% % get the location of the script file to save figures
% full_path = mfilename('fullpath');
% [startpath,  filename, ext] = fileparts(full_path);
% 
% %% get the folder listing
% 
% startpath = 'D:\IUPUI\Test_Data\';
% 
% data_path = uigetdir(startpath,'Select Scenario Folder');
% 
% if(data_path == 0)
%     return;
% end
% commandwindow;
% % 
% % listing = dir(data_path);
% % listing = listing(3:end);

%% cycle through each folder and try to 

% folders = strsplit(data_path, filesep);
% scenario_name = folders{end};
% camera_side = {'left','right'};


function get_sharpest_image(data_path, folders, scenario_name, camera_side)

    % cycle through each camera side and get the sharpest image for each exposure
    for idx=1:numel(camera_side)

        file_id = fopen(fullfile(data_path,strcat(scenario_name,'_',camera_side{idx},'_focus_file_listing.txt')),'w');

        fprintf('%s/\n',fullfile(folders{1:end-1}));
        fprintf(file_id,'# Data file used to select the infocus image\n%s/\n',fullfile(folders{1:end-1}));

        fprintf('#---------------------------------------------------------\n');
        fprintf(file_id,'#---------------------------------------------------------\n');    

        exp_listing = dir(strcat(data_path,'/',camera_side{idx},'/exp*'));
        %folders = strsplit(exp_listing(1).folder,filesep);

        for jdx=1:numel(exp_listing)

            exposure_name = exp_listing(jdx).name;
            full_data_path = strcat(scenario_name,'/',camera_side{idx},'/',exposure_name);
            image_listing = dir(strcat(fullfile(exp_listing(jdx).folder,exp_listing(jdx).name),'\*.png'));
            [img_vals] = get_sum_fft(image_listing);

            [max_val,index] = max(img_vals);
            %image_name = fullfile(image_listing(index).folder,image_listing(index).name);
            fprintf('%s/%s\n', full_data_path, image_listing(index).name);
            fprintf(file_id,'%s/%s\n', full_data_path, image_listing(index).name);
            fprintf(file_id,'#');
            for mdx=1:numel(img_vals)
                fprintf(file_id,'%f,',img_vals(mdx));
            end
            fprintf(file_id,'\n');            

        end

        fprintf('#---------------------------------------------------------\n');
        fprintf(file_id,'#---------------------------------------------------------\n');

        fclose(file_id);    

    end

end
