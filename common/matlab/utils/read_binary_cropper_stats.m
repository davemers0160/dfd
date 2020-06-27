function [data] = read_binary_cropper_stats(input_file)
    
    % open the file and read as binary
    fileID = fopen(input_file, 'r');
    
    % set up the waitbar
%     f = waitbar(0,'1','Name','',...
%         'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%     setappdata(f,'canceling',0);
    
    % get the file size
    fseek(fileID, 0, 'eof');
    file_size = ftell(fileID);
    fseek(fileID, 0, 'bof');
    
%     position = ftell(fileID);
    
    r = floor(file_size/24);
    
%     tmp = fread(fileID, file_size, 'uint32');
    
    %tic;
    tmp = fread(fileID, [6,r], 'uint32')';
    %toc;
    
    data = tmp(:,[1,3:6]);
    
%     index=1;
%     while (position < file_size)
%         
%         % Update waitbar and message
%         waitbar(position/file_size, f, 'Loading Data...');
%         
%         img_index = fread(fileID, 1, 'uint64');
%         tmp = fread(fileID, 4, 'uint32');
%         
%         if(~isempty(img_index) && numel(tmp)==4)
%             data(index,:) = [img_index, tmp'];
%             position = ftell(fileID);
%             index = index + 1;
%         else
%             fseek(fileID, 0,'eof');
%             break;
%         end
%     end
    
    fclose(fileID);
    
    %pause(0.2);
    %delete(f);
end