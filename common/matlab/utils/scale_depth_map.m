function scale_depth_map(listing, scale, save_path)

    for idx=1:length(listing)
        
        img_file = fullfile(listing(idx).folder, listing(idx).name);
        fprintf('%s\n', img_file);
        
        img = imread(img_file);
        
        img = scale*img;
        
        imwrite(uint8(img), fullfile(save_path, listing(idx).name));
        
    end

end
