% format long g
% format compact
% clc
% close all
% clearvars
% 
% full_path = mfilename('fullpath');
% [startpath,  filename, ext] = fileparts(full_path);
% plot_num = 1;

% %% Data Processing Main Code
% % This code will ask the user to select a config file that was generated
% % by the GUI to quickly parse all of the images in the raw data
% 
% % select a config file - the name of the file will dictate which scenario
% % and camera side will be parsed
% 
% file_filter = {'*.txt','Text Files';'*.*','All Files' };
% 
% startpath = 'D:\IUPUI\Test_Data\real_world_raw\configs';
% [cfg_file, cfg_file_path] = uigetfile(file_filter, 'Select Configuration File', startpath);
% if(cfg_file_path == 0)
%     return;
% end
% 
% % get the path to save the results
% startpath = 'D:\IUPUI\Test_Data\rw';
% save_path = uigetdir(startpath,'Select Save Folder');
% 
% if(save_path == 0)
%     return;
% end
% 
% load('d:\Common\matlab\lidar\OS1_991827000195.mat');
%% parse the name of the file
function process_data(cfg_file, cfg_file_path, v_step, save_path)

    cfg_parts = strsplit(cfg_file,{'_','.'});

    scenario_name = cfg_parts{2};
    camera_side = cfg_parts{3};
    
    folders = strsplit(cfg_file_path,filesep);
    data_path = fullfile(folders{1:end-1},scenario_name);

    %% parse through the
    focus_params = parse_input_parameters(fullfile(data_path, strcat(scenario_name,'_',camera_side,'_focus_file_listing.txt')));
    
    focus_params(1) = [];
    
    %% parse the contents of the file

    [match_params] = parse_input_parameters(fullfile(cfg_file_path, cfg_file));

    img_x_off = str2double(match_params{1}{1});
    img_y_off = str2double(match_params{1}{2});
    img_w = str2double(match_params{1}{3});
    img_h = str2double(match_params{1}{4});
    img_rot = str2double(match_params{1}{5});

    img_crop_w = [img_x_off:img_x_off+img_w-1];
    img_crop_h = [img_y_off:img_y_off+img_h-1];
    
    
%     ld_x_off = str2double(match_params{2}{1});
%     ld_y_off = str2double(match_params{2}{2});
%     ld_w = str2double(match_params{2}{3});
%     ld_h = str2double(match_params{2}{4});
% 
%     ld_tlx = str2double(match_params{3}{1});
%     ld_tly = str2double(match_params{3}{2});
%     ld_trx = str2double(match_params{3}{3});
%     ld_try = str2double(match_params{3}{4});
%     ld_blx = str2double(match_params{3}{5});
%     ld_bly = str2double(match_params{3}{6});
%     ld_brx = str2double(match_params{3}{7});
%     ld_bry = str2double(match_params{3}{8});
% 
%     ld_crop_w = [ld_x_off:ld_x_off+ld_w-1];
%     ld_crop_h = [ld_y_off:ld_y_off+ld_h-1];

    %% cycle through the exposure directories

%     save_focus_file_name = fullfile(data_path,strcat(scenario_name,'_',camera_side,'_input.txt'));
%     fileid = fopen(save_focus_file_name, 'w');
    
    exp_listing = dir(strcat(data_path,filesep,camera_side,filesep,'exp*'));

%     v_step_down = 5;
%     v_step_up = 0;
%     v_step_f = 133;
    
    for idx=1:length(exp_listing)

        exp_name = exp_listing(idx).name;
        
        tmp = strsplit(focus_params{idx}{1},'/');
        focus_file = tmp{end};
        
        % get the voltage step
        img_parts = strsplit(focus_file,'_');
        %v_step = str2double(img_parts{2});
        
        %v_step_min = min(v_step - v_step_down, v_step_f - v_step_down);
        %v_step_max = max(v_step + v_step_up, v_step_f + v_step_up);
        v_step_min = v_step(1);
        v_step_max = v_step(2);
        
        % get image listing
        tmp_listing = dir(strcat(exp_listing(idx).folder, filesep, exp_name, filesep,'*.png'));
        kdx = 1;
        clear img_listing;
        for jdx=1:length(tmp_listing)
            img_parts = strsplit(tmp_listing(jdx).name, {'_'});
            current_step = str2double(img_parts{2});
            
            if((current_step >= v_step_min) && (current_step <= v_step_max))
                img_listing(kdx,1) = tmp_listing(jdx);
                kdx = kdx + 1;
            end
        end
        
        for jdx=1:length(img_listing)

            img_parts = strsplit(img_listing(jdx).name, {'_'});
            img_name = strcat(img_parts{1}, '_', img_parts{2}, '_', img_parts{3}, '.png');

            img = imread(fullfile(img_listing(jdx).folder, img_listing(jdx).name));

            if(img_rot ~= 0.0)
                img_r = imrotate(img, img_rot, 'bilinear');
                img_cr = img_r(img_crop_h, img_crop_w, :);
            else
                img_cr = img(img_crop_h, img_crop_w, :);
            end

            save_location = fullfile(save_path, scenario_name, camera_side, exp_name);
            warning('off','all')
            mkdir(save_location);
            fprintf('Saving: %s\n', fullfile(save_location, img_name));
            imwrite(img_cr, fullfile(save_location, img_name));
            warning('on','all')
        end

    end

    %% save the lidar data 
    %process_lidar_data(lidar_struct, match_params, scenario_name, camera_side, data_path, save_path);
    
    % this is the scale factor to reduce the lidar data from a max of
%     % 80000mm to a number that fits into N-Bits starting at 8-bits
%     scale_factor = [320;160;80;40;20;10;5;2.5];
% 
%     save_location = fullfile(save_path,scenario_name,'lidar');
%     warning('off','all')
%     mkdir(save_location);
%     warning('on','all')
% 
%     ld_listing = dir(strcat(data_path, '/lidar/lidar_rng_00000*'));
%     lidar_file = ld_listing.name;
% 
%     ld_parts = strsplit(lidar_file,'_');
%     lidar_save_file = strcat(ld_parts{1},'_',ld_parts{2},'_',camera_side,'_',ld_parts{3},'.bin');
%     
%     lidar_save_file_png = {};
%     for idx=8:1:15
%         lidar_save_file_png{end+1} = strcat(ld_parts{1},'_',ld_parts{2},'_',camera_side,'_',ld_parts{3},'_', num2str(idx),'bit.png');
%     end
%     
%     %lidar_save_file_16_png = strcat(ld_parts{1},'_',ld_parts{2},'_',camera_side,'_',ld_parts{3},'_16bit.png');
%     %lidar_save_file_8_png = strcat(ld_parts{1},'_',ld_parts{2},'_',camera_side,'_',ld_parts{3},'_8bit.png');
% 
%     [ld_data] = read_binary_lidar_data(fullfile(ld_listing.folder, ld_listing.name));
%     [lx, ly, lz] = convert_lidar_to_xyz(lidar_struct, ld_data);
%     lx = lx(:,850:1150);
% 
%     [lx_height, lx_width] = size(lx);
%     max_data = (floor(max(lx(:))/1000) + 1)*1000;
%     mp = [1,1; lx_width,1; 1, lx_height; lx_width, lx_height];
%     fp = [mp(1,1)+ld_tlx, mp(1,2)+ld_tly; mp(2,1)+ld_trx, mp(2,2)+ld_try; mp(3,1)+ld_blx, mp(3,2)+ld_bly; mp(4,1)+ld_brx, mp(4,2)+ld_bry];
% 
%     tform = fitgeotrans(mp,fp,'projective');
%     RA = imref2d([lx_height lx_width], [1 lx_width], [1 lx_height]);
% 
%     [lx_w, ref] = imwarp(lx, tform, 'OutputView', RA);
%     ld_cr = max(0, lx_w(ld_crop_h, ld_crop_w));
%     ld_cr_msk = (ld_cr < 500); % find the values that are less than 500 mm
% 
% 
%     fprintf('Saving: %s\n', fullfile(save_location, lidar_save_file));
%     write_uint32_binary_image(fullfile(save_location, lidar_save_file), ld_cr);
%     
%     for idx = 1:numel(lidar_save_file_png)
%         %ld_cr_bit = uint16(floor(ld_cr/(2^(20-(7+idx)))));     % gives a true bit reduction based on a 20-bit number
%         ld_cr_bit = uint16(floor(ld_cr/scale_factor(idx)));          % reduces the number based on a max of 80000mm
%         ld_cr_bit(ld_cr_msk) = 65535;
%         fprintf('Saving: %s\n', fullfile(save_location, lidar_save_file_png{idx}));
%         imwrite(ld_cr_bit, fullfile(save_location, lidar_save_file_png{idx}), 'BitDepth', 16);
%     end

%     fprintf('Saving: %s\n', fullfile(save_location, lidar_save_file_16_png));
%     imwrite(ld_cr_16, fullfile(save_location, lidar_save_file_16_png), 'BitDepth', 16);
% 
%     fprintf('Saving: %s\n', fullfile(save_location, lidar_save_file_8_png));
%     imwrite(ld_cr_8, fullfile(save_location, lidar_save_file_8_png), 'BitDepth', 8);

end
%fprintf('Complete!\n');

