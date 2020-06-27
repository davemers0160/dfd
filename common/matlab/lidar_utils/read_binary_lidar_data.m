
function [data] = read_binary_lidar_data(input_file)

    % open the file and read as a little-endian 
    file_id = fopen(input_file, 'r', 'l');
    
    % read the size
    sz = fread(file_id, 2, 'uint32');
     
    % read in the data
    data = fread(file_id, sz(1)*sz(2), 'int32');
    data = reshape(data,[sz(2),sz(1)])';

    data = single(data);
    fclose(file_id);
    
end