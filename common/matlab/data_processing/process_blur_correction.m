function process_blur_correction(cfg_file, cfg_file_path, v_step, save_path, kernel, mask_params)

    cfg_parts = strsplit(cfg_file,{'_','.'});

    scenario_name = cfg_parts{2};
    camera_side = cfg_parts{3};

    folders = strsplit(cfg_file_path, filesep);
    data_path = fullfile(folders{1:end-1}, scenario_name);

    %% parse the contents of the file

    [match_params] = parse_input_parameters(fullfile(cfg_file_path, cfg_file));

    img_x_off = str2double(match_params{1}{1});
    img_y_off = str2double(match_params{1}{2});
    img_w = str2double(match_params{1}{3});
    img_h = str2double(match_params{1}{4});
    img_rot = str2double(match_params{1}{5});

    img_crop_w = [img_x_off:img_x_off+img_w-1];
    img_crop_h = [img_y_off:img_y_off+img_h-1];

    % get the exposre directory listing
    exp_listing = dir(strcat(data_path,filesep,camera_side,filesep,'exp*'));

    for idx=1:numel(exp_listing)

        exp_name = exp_listing(idx).name;
        img_listing=[];
        img_listing = dir(strcat(exp_listing(idx).folder, filesep, exp_name, filesep,'*',num2str(v_step(1)),'*.png'));

        img_parts = strsplit(img_listing.name, {'_'});

        defocus_img = double(imread(fullfile(img_listing(1).folder,img_listing(1).name)));

        [corr_img] = apply_blur_correction(defocus_img, kernel, mask_params);

        save_img_name = strcat(img_parts{1}, '_', img_parts{2}, '_', img_parts{3}, '_corr.png');

        if(img_rot ~= 0.0)
            img_r = imrotate(corr_img, img_rot, 'bilinear');
            img_cr = img_r(img_crop_h, img_crop_w, :);
        else
            img_cr = corr_img(img_crop_h, img_crop_w, :);
        end

        save_location = fullfile(save_path, scenario_name, camera_side, exp_name);
        warning('off','all')
        mkdir(save_location);
        fprintf('Saving: %s\n', fullfile(save_location, save_img_name));
        imwrite(uint8(img_cr), fullfile(save_location, save_img_name));
        warning('on','all')
    end

end