
function [log_version, data, final_traning, final_test] = read_dnn_logfile(log_file)

    sep = '------------------------------------------------------------------';

    file_id = fopen(log_file,'r');

    index = 1;
    tline = fgetl(file_id); % read in the first line
    tline = fgetl(file_id); % read in the second line and get the version number

    l = strsplit(tline,{':','.',' '});

    log_version.major = str2double(l{2});
    log_version.minor = str2double(l{3});

    %% run through the logs until we get to the results section
    while ~feof(file_id)

        tline = fgetl(file_id);

        if(strcmp(tline,sep))
            tline = fgetl(file_id);
            l = strsplit(tline,{','});

            if(strcmp(l{1},'<step>'))
                break;            
            end
        end    
    end

    %% start getting the data
    %<step>, <learning rate>, <average loss>, Train (<SAE>, <MSE>, <SSIM>), Test (<SAE>, <MSE>, <SSIM>)
    data = [];

    while ~feof(file_id)
        tline = fgetl(file_id);
        l = strsplit(tline,{',',':'});

        if(~strcmp(l{1},'Stop Code'))
            data(index).step = str2double(l{1});
            data(index).lr = str2double(l{2});
            data(index).avg_loss = str2double(l{3});
            data(index).tr_nmae = str2double(l{4});
            data(index).tr_nrmse = str2double(l{5});
            data(index).tr_ssim = str2double(l{6});
            data(index).te_nmae = str2double(l{7});
            data(index).te_nrmse = str2double(l{8});
            data(index).te_ssim = str2double(l{9});

            index = index + 1;
        else
            break;
        end
    end

    %% get the final training results
    final_traning = [];

    while ~feof(file_id)
        tline = fgetl(file_id);

        if(strcmp(tline,sep))
            tline = fgetl(file_id);
            tline = fgetl(file_id);
            l = strsplit(tline,{',',':'});

            if(numel(l)==6)
                final_traning(1) = str2double(l{4});
                final_traning(2) = str2double(l{5});
                final_traning(3) = str2double(l{6});
                break;
            end

        end  
    end

    %% get the final testing results
    final_test = [];

    while ~feof(file_id)
        tline = fgetl(file_id);

        if(strcmp(tline,sep))
            tline = fgetl(file_id);
            tline = fgetl(file_id);
            l = strsplit(tline,{',',':'});

            if(numel(l)==6)
                final_test(1) = str2double(l{4});
                final_test(2) = str2double(l{5});
                final_test(3) = str2double(l{6});
                break;
            end

        end  
    end
    
    %% close the file
    fclose(file_id);

end

