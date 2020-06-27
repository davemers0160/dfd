% get the file
file_filter = {'*.txt';'*.*'};
save_path = 'D:\IUPUI\DfD\DFD_GC';
[input_file, save_path] = uiputfile(file_filter, 'Enter save file', save_path);

if(save_path == 0)
    return;
end

%%
data_directory = 'D:\IUPUI\Test_Data\rw';
scenario = 'Auditorium2';
exp = {'exp_10', 'exp_20', 'exp_30', 'exp_40', 'exp_50', 'exp_60', 'exp_70'};

min_step = 128;
max_step = 142;
v=(min_step:1:max_step);
nk = nchoosek(v,2);

file_id = fopen(fullfile(save_path,input_file),'w');

fprintf(file_id, '%s\n', data_directory);
fprintf(file_id, '#data: \n');

for idx=1:size(nk,1)

    fprintf('%s/left/%s/image_%d_50.00.png, %s/left/%s/image_%d_50.00.png, %s/lidar/lidar_rng_left_00000_8bit.png\n', ...
        scenario, exp{5}, nk(idx,1), scenario, exp{5}, nk(idx,2), scenario);    
    fprintf(file_id,'%s/left/%s/image_%d_50.00.png, %s/left/%s/image_%d_50.00.png, %s/lidar/lidar_rng_left_00000_8bit.png\n', ...
        scenario, exp{5}, nk(idx,1), scenario, exp{5}, nk(idx,2), scenario);
%     fprintf(file_id,'%s/right/%s/image_%d_50.00.png, %s/right/%s/image_%d_50.00.png, %s/lidar/lidar_rng_right_00000_8bit.png\n', ...
%         scenario, exp{5}, nk(idx,1), scenario, exp{5}, nk(idx,2), scenario);
end

fclose(file_id);