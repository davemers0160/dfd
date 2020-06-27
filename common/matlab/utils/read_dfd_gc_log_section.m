function [data, data_name, image_files] = read_dfd_gc_log_section(file_id, sep)
    data_name = {};
    image_files = {};
    data = [];
    
    tline = fgetl(file_id); % read in the second line - focus file name
    while(isempty(tline) && ~feof(file_id))
        tline = fgetl(file_id); % read in the second line - focus file name     
    end
    
    if(feof(file_id))
        return;
    end
    
    % parse the path into its folder parts
    folders = regexp(tline,'/','split');
    
    % the structure will be in two parts
    % the name will be a cell array of string values
    % the data will be in a numerical array
    data_name = folders{numel(folders)-3};
    
    infocus_img_file = fullfile(folders{numel(folders)-3:numel(folders)});

    switch folders{numel(folders)-2}
        case 'Illum1'
            I = 1;
        case 'Illum2'
            I = 2;
        case 'Illum3'
            I = 3;
    end

    switch folders{numel(folders)-1}
        case 'Exp0'
            E = 1;
        case 'Exp1'
            E = 2;
        case 'Exp2'
            E = 3; 
    end

    switch folders{numel(folders)}
        case 'view1.png'
            V = 1;
        case 'view5.png'
            V = 2;
    end

    % get the defocused image location
    tline = fgetl(file_id);
    folders = regexp(tline,'/','split');
    defocus_img_file = fullfile(folders{numel(folders)-3:numel(folders)});

    % get the defocused image location
    tline = fgetl(file_id);
    folders = regexp(tline,'/','split');
    gt_img_file = fullfile(folders{numel(folders)-1:numel(folders)});

    run_time = 0;
    while ~feof(file_id)

        tline = fgetl(file_id);

        s_line = regexp(tline,': ','split');

        if(strcmp(s_line{1}, 'Depth Map Generation Complete (minutes)'))
            run_time = str2double(s_line{2});
        end

        if(strcmp(tline,sep))
            break;
        end
    end

    tline = fgetl(file_id);

    results = regexp(tline,'\t','split');
    r2 = regexp(results{1},': ','split');
    r2_2 = regexp(r2{2},', ', 'split');
    
    NMAE = str2double(r2_2{1});
    NRMSE = str2double(r2_2{2});
    SSIM = str2double(r2_2{3});

    data = [I, E, V, NMAE, NRMSE, SSIM, run_time];
    
    image_files = {infocus_img_file, defocus_img_file, gt_img_file};
    
    % read the last separator line
    tline = fgetl(file_id);

end