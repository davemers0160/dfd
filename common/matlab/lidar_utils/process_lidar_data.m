function process_lidar_data(lidar_struct, match_params, scenario_name, camera_side, data_path, save_path)

    % 80000mm to a number that fits into N-Bits starting at 8-bits
    %scale_factor = [320;160;80;40;20;10;5;2.5];
    scale_factor = [10;5;4;2;1;1;1;1];

    ld_x_off = str2double(match_params{2}{1});
    ld_y_off = str2double(match_params{2}{2});
    ld_w = str2double(match_params{2}{3});
    ld_h = str2double(match_params{2}{4});

    ld_tlx = str2double(match_params{3}{1});
    ld_tly = str2double(match_params{3}{2});
    ld_trx = str2double(match_params{3}{3});
    ld_try = str2double(match_params{3}{4});
    ld_blx = str2double(match_params{3}{5});
    ld_bly = str2double(match_params{3}{6});
    ld_brx = str2double(match_params{3}{7});
    ld_bry = str2double(match_params{3}{8});
      
    ld_offset = zeros(ld_h,1);
    if(size(match_params,1)>=4)       
        for idx=1:numel(match_params{4})
            ld_offset(idx,1) = str2double(match_params{4}{idx});
        end
    end

    ld_crop_w = [ld_x_off:ld_x_off+ld_w-1];
    ld_crop_h = [ld_y_off:ld_y_off+ld_h-1];
    
    save_location = fullfile(save_path,scenario_name,'lidar');
    warning('off','all')
    mkdir(save_location);
    warning('on','all')

    ld_listing = dir(strcat(data_path, '/lidar/lidar_rng_00000*'));
    lidar_file = ld_listing.name;

    ld_parts = strsplit(lidar_file,'_');
    lidar_save_file = strcat(ld_parts{1},'_',ld_parts{2},'_',camera_side,'_',ld_parts{3},'.bin');
    
    lidar_save_file_png = {};
    for idx=8:1:15
        lidar_save_file_png{end+1} = strcat(ld_parts{1},'_',ld_parts{2},'_',camera_side,'_',ld_parts{3},'_', num2str(idx),'bit.png');
    end
    
    [ld_data] = read_binary_lidar_data(fullfile(ld_listing.folder, ld_listing.name));
    [lx, ly, lz] = convert_lidar_to_xyz(lidar_struct, ld_data);
    lx = lx(:,850:1150);

    [lx_height, lx_width] = size(lx);
    %max_data = (floor(max(lx(:))/1000) + 1)*1000;
    
    mp = [1,1; lx_width,1; 1, lx_height; lx_width, lx_height];
    fp = [mp(1,1)+ld_tlx, mp(1,2)+ld_tly; mp(2,1)+ld_trx, mp(2,2)+ld_try; mp(3,1)+ld_blx, mp(3,2)+ld_bly; mp(4,1)+ld_brx, mp(4,2)+ld_bry];

    tform = fitgeotrans(mp,fp,'projective');
    RA = imref2d([lx_height lx_width], [1 lx_width], [1 lx_height]);

    [lx_w, ref] = imwarp(lx, tform, 'OutputView', RA);
    %ld_cr = max(0, lx_w(ld_crop_h, ld_crop_w));
    
    if(lx_w == 0)
        ld_cr = max(0, lx_w(ld_crop_h, ld_crop_w));
    else  
        for r=1:ld_h   
           if(ld_offset(r)>0)
               tmp_lx(r,:) = lx_w(r + ld_y_off - 1, ...
                   [lx_width-ld_offset(r):lx_width, 1:(lx_width-ld_offset(r)-1)]);
           elseif(ld_offset(r) < 0)
               tmp_lx(r,:) = lx_w(r + ld_y_off - 1, ...
                   [abs(ld_offset(r)):lx_width, 1:abs(ld_offset(r))-1]);
           else
               tmp_lx(r,:) = lx_w(r + ld_y_off - 1,:);
           end
        end    
        ld_cr = tmp_lx(:, ld_crop_w);
    end

    %ld_cr_msk = (ld_cr < 500); % find the values that are less than 500 mm

    fprintf('Saving: %s\n', fullfile(save_location, lidar_save_file));
    write_uint32_binary_image(fullfile(save_location, lidar_save_file), ld_cr);
    
    for idx = 1:numel(lidar_save_file_png)
        %ld_cr_bit = uint16(floor(ld_cr/(2^(20-(7+idx)))));     % gives a true bit reduction based on a 20-bit number
        ld_cr_bit = uint16(floor(ld_cr/scale_factor(idx)));          % reduces the number based on a max of 80000mm
        %ld_cr_bit(ld_cr_msk) = 65535;
        fprintf('Saving: %s\n', fullfile(save_location, lidar_save_file_png{idx}));
        imwrite(ld_cr_bit, fullfile(save_location, lidar_save_file_png{idx}), 'BitDepth', 16);
    end     
    
end
