function [data] = read_binary_depthmap_file(input_file, size_type, data_type)

    % open the file and read as a little-endian 
    file_id = fopen(input_file, 'r', 'l');
    
    % read the size
    sz = fread(file_id, 2, size_type);
     
    % read in the data
    data = fread(file_id, sz(1)*sz(2), data_type);
    data = reshape(data,[sz(2),sz(1)])';

    %data = single(data);
    fclose(file_id);
    
    
end