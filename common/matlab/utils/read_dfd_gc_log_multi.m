function [data, data_name, image_files] = read_dfd_gc_log_multi(logfile, listing, sep)

    data_name = {};
    image_files = {};
    data = [];

    if(isempty(logfile))
       fprintf('No log file found...\n');
       return;
    end

    dm_file = dir(strcat(fullfile(listing.folder,listing.name),'\*.png'));
    if (isempty(dm_file))
        dm_img_file = {''};
    else
        dm_img_file = {fullfile(dm_file(2).folder,dm_file(2).name)};
    end
    
    file_id = fopen(fullfile(logfile(1).folder,logfile(1).name),'r');

    platform = fgetl(file_id); % read in the fisrt line - platform
    
    
    while(~feof(file_id))
      
        [d, d_n, i_f] = read_dfd_gc_log_section(file_id, sep);

        i_f = cat(2, i_f, dm_img_file);
    
        if(~isempty(d))
            data(end+1,:) = d;
            data_name{end+1,1} = d_n;
            image_files{end+1,1} = i_f; 
        end
    end

    fclose(file_id);

end