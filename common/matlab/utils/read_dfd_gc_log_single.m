function [data, data_name, image_files] = read_dfd_gc_log_single(logfile, listing, sep)

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
      
    [data, data_name, image_files] = read_dfd_gc_log_section(file_id, sep);
    
    image_files = cat(2, image_files, dm_img_file);

    fclose(file_id);

end