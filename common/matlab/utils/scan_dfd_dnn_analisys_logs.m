format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

% This scripts now only reads the  analysis log files v2.1+

%% setup some variables
plot_num = 1;
sep = '------------------------------------------------------------------';
file_filter = {'*.txt','Text Files';'*.*','All Files' };

generate_images = false;

%% get the directory to scan
startpath = 'D:\IUPUI\PhD\Results';
[log_file, log_path] = uigetfile(file_filter, 'Select DFD-Net Log File', startpath);
if(log_path == 0)
    return;
end

commandwindow;

%% Start reading the file
data_name = {};
focus_file = {};
defocus_file = {};
data = [];

file_id = fopen(fullfile(log_path,log_file),'r');

index = 1;

while ~feof(file_id)
    
    tline = fgetl(file_id);
    
    if(strcmp(tline,sep))
        tline = fgetl(file_id);
        
        if(~isempty(tline) && strcmp(tline(1:8),'Depthmap'))
            
            % run time
            r_line = regexp(tline,': ','split');
            r_line2 = regexp(r_line{2},' ','split');
            run_time = str2double(r_line2{1});
            
            % image size
            tline = fgetl(file_id);
            i_line = regexp(tline,': ','split');
            i_line = regexp(i_line{2},' x ', 'split');
            img_size(index,:) = [str2double(i_line{1}), str2double(i_line{2})];
            
            % focus file name
            tline = fgetl(file_id);
            f_line = regexp(tline,'Focus File: ', 'split');
            folders = regexp(f_line{2},'/','split');
            focus_file{index,1} = fullfile(folders{end-3:end});
            
            % ground truth file 
            gt_file{index,1} = fullfile(folders{end-3},strcat('disp',folders{end}(5),'.png'));
            
            data_name{index,1} = folders{numel(folders)-3};
            
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
            
            % defocus file
            tline = fgetl(file_id);
            f_line = regexp(tline,'Defocus File: ','split');
            folders = regexp(f_line{2},'/','split');
            defocus_file{index,1} = fullfile(folders{end-3:end});  
            
            % depth map file
            tline = fgetl(file_id);
            f_line = regexp(tline,': ','split');
            folders = regexp(f_line{2},'/','split');
            dm_file{index,1} = fullfile(folders{end}); 
            
            tline = fgetl(file_id);
            r2 = regexp(tline,': ','split');
            NMAE = str2double(r2{2});
            
            tline = fgetl(file_id);
            r2 = regexp(tline,': ','split');
            NRMSE = str2double(r2{2});
            
            tline = fgetl(file_id);
            r2 = regexp(tline,': ','split');
            SSIM = str2double(r2{2});
            
%             tline = fgetl(file_id);
%             r2 = regexp(tline,': ','split');
%             SILOG = str2double(r2{2});            
            
%            data(index,:) = [I, E, V, NMAE, NRMSE, SSIM, SILOG, run_time];
            data(index,:) = [I, E, V, NMAE, NRMSE, SSIM, 0, run_time];
            index = index + 1;
            
        end

    end    
    
end

clear NRMSE NMAE SSIM SILOG

%% close the file

fclose(file_id);

%%
%indexes for the data
iNMAE = 4;
iNRMSE = 5;
iSSIM = 6;
iSILOG = 7;
irun = 8;
scenario = 'test';

% do this to get a baseline set of names for the xticklabel and to fill in
% the blanks for datasets with no data
plot_I = data(data(:,1) == 1,:);
I_name = data_name(data(:,1) == 1,:);
x_name = I_name(plot_I(:,2) == 1,:);
cm = hsv(3);

exposure_legend = {strcat('Exposure 0',32), strcat('Exposure 1',32), strcat('Exposure 2',32)};

for illum=1:3
    plot_I = data(data(:,1) == illum,:);
    I_name = data_name(data(:,1) == illum,:);

    plot_E = {};
    E_name = {};
    
    for idx=1:3
        plot_E{idx} = plot_I(plot_I(:,2) == idx,:);
        E_name{idx} = I_name(plot_I(:,2) == idx,:);
    end
    
    avg_runtime(illum,1) = mean(plot_E{1}(:,irun));
    avg_runtime(illum,2) = mean(plot_E{2}(:,irun));
    avg_runtime(illum,3) = mean(plot_E{3}(:,irun));   
    
%     NRMSE_mean(illum,1) = mean(plot_E{1}(:,iNRMSE),'omitnan');
%     NRMSE_mean(illum,2) = mean(plot_E{2}(:,iNRMSE),'omitnan');
%     NRMSE_mean(illum,3) = mean(plot_E{3}(:,iNRMSE),'omitnan');
%     
%     NRMSE_std(illum,1) = std(plot_E{1}(:,iNRMSE),'omitnan');
%     NRMSE_std(illum,2) = std(plot_E{2}(:,iNRMSE),'omitnan');
%     NRMSE_std(illum,3) = std(plot_E{3}(:,iNRMSE),'omitnan');
% 
%     NMAE_mean(illum,1) = mean(plot_E{1}(:,iNMAE),'omitnan');
%     NMAE_mean(illum,2) = mean(plot_E{2}(:,iNMAE),'omitnan');
%     NMAE_mean(illum,3) = mean(plot_E{3}(:,iNMAE),'omitnan');
%     
%     NMAE_std(illum,1) = std(plot_E{1}(:,iNMAE),'omitnan');
%     NMAE_std(illum,2) = std(plot_E{2}(:,iNMAE),'omitnan');
%     NMAE_std(illum,3) = std(plot_E{3}(:,iNMAE),'omitnan');
%     
%     SSIM_mean(illum,1) = mean(plot_E{1}(:,iSSIM),'omitnan');
%     SSIM_mean(illum,2) = mean(plot_E{2}(:,iSSIM),'omitnan');
%     SSIM_mean(illum,3) = mean(plot_E{3}(:,iSSIM),'omitnan');
%     
%     SSIM_std(illum,1) = std(plot_E{1}(:,iSSIM),'omitnan');
%     SSIM_std(illum,2) = std(plot_E{2}(:,iSSIM),'omitnan');
%     SSIM_std(illum,3) = std(plot_E{3}(:,iSSIM),'omitnan');    
    
    for idx=1:numel(x_name)
        if(strcmp(x_name{idx},E_name{1,1}{idx})==0)
            E_name{1,1} = [{E_name{1,1}{1:idx-1}}';{x_name{idx}};{E_name{1,1}{idx:end}}'];
            plot_E{1,1} = [plot_E{1,1}(1:idx-1,:); [illum, 2, nan, nan, nan, nan, nan, nan]; plot_E{1,1}(idx:end,:)];
        end         
        if(strcmp(x_name{idx},E_name{1,2}{idx})==0)
            E_name{1,2} = [{E_name{1,2}{1:idx-1}}';{x_name{idx}};{E_name{1,2}{idx:end}}'];
            plot_E{1,2} = [plot_E{1,2}(1:idx-1,:); [illum, 2, nan, nan, nan, nan, nan, nan]; plot_E{1,2}(idx:end,:)];
        end    
        if(strcmp(x_name{idx},E_name{1,3}{idx})==0)
            E_name{1,3} = [{E_name{1,3}{1:idx-1}}';{x_name{idx}};{E_name{1,3}{idx:end}}'];
            plot_E{1,3} = [plot_E{1,3}(1:idx-1,:); [illum, 3, nan, nan, nan, nan, nan, nan]; plot_E{1,3}(idx:end,:)];
        end
    end   

    NRMSE(illum,:,1) = plot_E{1}(:,iNRMSE);
    NRMSE(illum,:,2) = plot_E{2}(:,iNRMSE);
    NRMSE(illum,:,3) = plot_E{3}(:,iNRMSE);
    
    NMAE(illum,:,1) = plot_E{1}(:,iNMAE);
    NMAE(illum,:,2) = plot_E{2}(:,iNMAE);
    NMAE(illum,:,3) = plot_E{3}(:,iNMAE);

    SSIM(illum,:,1) = plot_E{1}(:,iSSIM);
    SSIM(illum,:,2) = plot_E{2}(:,iSSIM);
    SSIM(illum,:,3) = plot_E{3}(:,iSSIM);

    SILOG(illum,:,1) = plot_E{1}(:,iSILOG);
    SILOG(illum,:,2) = plot_E{2}(:,iSILOG);
    SILOG(illum,:,3) = plot_E{3}(:,iSILOG);
    
%     % plot the NRMSE and SSIM
%     nrmse_max = max(data(:,iNRMSE));
%     fig = figure(plot_num);
%     set(gcf,'position',([100,100,1000,600]),'color','w')
%     set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
%     hold on
%     box on
%     grid on
% 
%     set(gca,'fontweight','bold')
%     yyaxis left
%     plot(plot_E{1}(:,iNRMSE))
%     plot(plot_E{2}(:,iNRMSE))
%     plot(plot_E{3}(:,iNRMSE))
%     ylim([0 ceil(nrmse_max)])
%     ylabel('NRMSE','fontweight','bold')
%     %title(strcat(10,'Illumination',32,num2str(illum),10))
% 
%     yyaxis right
%     plot(plot_E{1}(:,iSSIM))
%     plot(plot_E{2}(:,iSSIM))
%     plot(plot_E{3}(:,iSSIM))
%     ylim([0.0 1.0])
%     ylabel('SSIM','fontweight','bold')
% 
%     xlim([1 length(E_name{1})]);
%     xticks([1:length(E_name{1})]);
%     xticklabels(E_name{1});
%     xtickangle(90);
%     legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
%     ax = gca;
%     ax.Position = [0.07 0.23 0.85 0.75];
% 
%     print(plot_num, '-dpng', fullfile(log_path,strcat('DNN_NRMSE_Results_I',num2str(illum),'_',scenario,'.png')));
% 
%     plot_num = plot_num + 1;
% 
%     % plot the NMAE and SSIM
% 
%     mae_max = max(data(:,iNMAE));
%     fig = figure(plot_num);
%     set(gcf,'position',([100,100,1000,600]),'color','w')
%     set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
%     hold on
%     box on
%     grid on
% 
%     set(gca,'fontweight','bold')
%     yyaxis left
%     plot(plot_E{1}(:,iNMAE))
%     plot(plot_E{2}(:,iNMAE))
%     plot(plot_E{3}(:,iNMAE))
%     ylim([0 ceil(mae_max)])
%     ylabel('NMAE','fontweight','bold')
%     %title(strcat(10,'Illumination',32,num2str(illum),10))
% 
%     yyaxis right
%     plot(plot_E{1}(:,iSSIM))
%     plot(plot_E{2}(:,iSSIM))
%     plot(plot_E{3}(:,iSSIM))
%     ylim([0.0 1.0])
%     ylabel('SSIM','fontweight','bold')
% 
%     xlim([1 length(E_name{1})]);
%     xticks([1:length(E_name{1})]);
%     xticklabels(E_name{1});
%     xtickangle(90);
%     legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
%     ax = gca;
%     ax.Position = [0.07 0.23 0.85 0.75];
% 
%     print(plot_num, '-dpng', fullfile(log_path,strcat('DNN_NMAE_Results_I',num2str(illum),'_',scenario,'.png')));
% 
%     plot_num = plot_num + 1;
    
    run_times{illum,1} = plot_E{1}(:,irun);
    run_times{illum,2} = plot_E{2}(:,irun);
    run_times{illum,3} = plot_E{3}(:,irun);
    
%     % plot the run times for each exposure level
%     figure(plot_num);
%     set(gcf,'position',([100,100,1000,600]),'color','w')
%     hold on
%     box on
%     grid on
%     plot(plot_E{1}(:,irun),'-b')
%     plot(plot_E{2}(:,irun),'--b')
%     plot(plot_E{3}(:,irun),':b')
%     set(gca,'fontweight','bold')
%     xlim([1 length(E_name{1})]);
%     xticks([1:length(E_name{1})]);
%     xticklabels(E_name{1});
%     xtickangle(90);
% 
%     ylabel('Run Time (s)', 'fontweight','bold')
%     legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
%     ax = gca;
%     ax.Position = [0.07 0.23 0.85 0.75];
%  
%     print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_DNN_Runtime_I',num2str(illum),'_',scenario,'.png')));
%     plot_num = plot_num + 1;
    
    % now plot everything in its own subplot just to see how it looks
    fig = figure(plot_num);
    set(gcf,'position',([50,50,1400,700]),'color','w')
    nrmse_max = max(data(:,iNRMSE));
    
    subplot(3,1,1)
    hold on
    box on
    grid on

%     plot(plot_E{1}(:,iNRMSE),'-b')
%     plot(plot_E{2}(:,iNRMSE),'--b')
%     plot(plot_E{3}(:,iNRMSE),':b')
    b1=bar([plot_E{1}(:,iNRMSE),plot_E{2}(:,iNRMSE),plot_E{3}(:,iNRMSE)]);
    set(gca,'fontweight','bold','FontSize', 13)
    
    % X-Axis
    xlim([0 length(E_name{1})+1]);
    xticks([1:length(E_name{1})]);
    xticklabels([]);
    xtickangle(90);
    
    % Y-Axis
    ylim([0 ceil(nrmse_max*10)/10])
    ylabel('NRMSE','fontweight','bold','FontSize', 13)
    
    title(strcat('DFD-Net Performance Results - Illumination',32,num2str(illum)),'fontweight','bold','FontSize', 16)

    b1(1).FaceColor = cm(1,:);
    b1(2).FaceColor = cm(2,:);
    b1(3).FaceColor = cm(3,:);  
    ax = gca;
    ax.Position = [0.05 0.715 0.93 0.24];
    
    nmae_max = max(data(:,iNMAE));    
    subplot(3,1,2)
    hold on
    box on
    grid on

%     plot(plot_E{1}(:,iNMAE),'-b')
%     plot(plot_E{2}(:,iNMAE),'--b')
%     plot(plot_E{3}(:,iNMAE),':b')
    b2=bar([plot_E{1}(:,iNMAE),plot_E{2}(:,iNMAE),plot_E{3}(:,iNMAE)]);
    set(gca,'fontweight','bold','FontSize', 13)
    
    % X-Axis
    xlim([0 length(E_name{1})+1]);
    xticks([1:length(E_name{1})]);
    xticklabels([]);
    xtickangle(90);  
    
    % X-Axis
    ylim([0 ceil(nmae_max*10)/10])
    ylabel('NMAE','fontweight','bold','FontSize', 13)


    b2(1).FaceColor = cm(1,:);
    b2(2).FaceColor = cm(2,:);
    b2(3).FaceColor = cm(3,:);  
    ax = gca;
    ax.Position = [0.05 0.445 0.93 0.24];

    subplot(3,1,3)
    hold on
    box on
    grid on

%     plot(plot_E{1}(:,iSSIM),'-b')
%     plot(plot_E{2}(:,iSSIM),'--b')
%     plot(plot_E{3}(:,iSSIM),':b')
    b3=bar([plot_E{1}(:,iSSIM),plot_E{2}(:,iSSIM),plot_E{3}(:,iSSIM)]);
    set(gca,'fontweight','bold','FontSize', 13)
    
    % X-Axis
    xlim([0 length(E_name{1})+1]);
    xticks([1:length(E_name{1})]);
    xticklabels(E_name{1});
    xtickangle(90);
    
    % Y-Axis
    ylim([0 1]);
    ylabel('SSIM','fontweight','bold','FontSize', 13);
    yticks((0:0.25:1));

    b3(1).FaceColor = cm(1,:);
    b3(2).FaceColor = cm(2,:);
    b3(3).FaceColor = cm(3,:);  
    %legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
    legend(exposure_legend, 'location', 'southoutside', 'orientation', 'horizontal')
    ax = gca;
    ax.Position = [0.05 0.175 0.93 0.24];

    print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_DNN_Perf_I',num2str(illum),'_',scenario,'.png')));

    plot_num = plot_num + 1;      
    
end

%% generate a runtime graph for each of the illuminations in subplots

run_max = max(data(:,irun));

figure(plot_num);
set(gcf,'position',([50,50,1400,700]),'color','w')

illum = 1;
subplot(3,1,1)
hold on
box on
grid on
% plot(run_times{illum,1},'-b')
% plot(run_times{illum,2},'--b')
% plot(run_times{illum,3},':b')
b1 = bar([run_times{illum,1},run_times{illum,2},run_times{illum,3}]);
set(gca,'fontweight','bold','FontSize', 13)

% X-Axis
xlim([0 length(E_name{1})+1]);
xticks([1:length(E_name{1})]);
xticklabels([]);
xtickangle(90);

% Y-Axis
ylim([0,ceil(run_max*10)/10]);
ylabel('Run Time (s)', 'fontweight','bold','FontSize', 13)

b1(1).FaceColor = cm(1,:);
b1(2).FaceColor = cm(2,:);
b1(3).FaceColor = cm(3,:); 

title(strcat('DFD-Net Run Time - Illumination',32,num2str(illum)),'fontweight','bold','FontSize', 13)
%legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.05 0.715 0.93 0.22];

illum = 2;
subplot(3,1,2)
hold on
box on
grid on
% plot(run_times{illum,1},'-b')
% plot(run_times{illum,2},'--b')
% plot(run_times{illum,3},':b')
b2 = bar([run_times{illum,1},run_times{illum,2},run_times{illum,3}]);
set(gca,'fontweight','bold','FontSize', 13)

% X-Axis
xlim([0 length(E_name{1})+1]);
xticks([1:length(E_name{1})]);
xticklabels([]);
xtickangle(90);

% Y-Axis
ylim([0,ceil(run_max*10)/10]);
ylabel('Run Time (s)', 'fontweight','bold','FontSize', 13)

b2(1).FaceColor = cm(1,:);
b2(2).FaceColor = cm(2,:);
b2(3).FaceColor = cm(3,:); 

title(strcat('DFD-Net Run Time - Illumination',32,num2str(illum)))
ax = gca;
ax.Position = [0.05 0.445 0.93 0.22];

illum = 3;
subplot(3,1,3)
hold on
box on
grid on
% plot(run_times{illum,1},'-b')
% plot(run_times{illum,2},'--b')
% plot(run_times{illum,3},':b')
b3 = bar([run_times{illum,1},run_times{illum,2},run_times{illum,3}]);
set(gca,'fontweight','bold','FontSize', 13)

% X-Axis
xlim([0 length(E_name{1})+1]);
xticks([1:length(E_name{1})]);
xticklabels(E_name{1});
xtickangle(90);

% Y-Axis
ylim([0,ceil(run_max*10)/10]);
ylabel('Run Time (s)', 'fontweight','bold','FontSize', 13)

title(strcat('DFD-Net Run Time - Illumination',32,num2str(illum)),'fontweight','bold','FontSize', 13)

b3(1).FaceColor = cm(1,:);
b3(2).FaceColor = cm(2,:);
b3(3).FaceColor = cm(3,:); 
%legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
legend(exposure_legend, 'location', 'southoutside', 'orientation', 'horizontal');

ax = gca;
ax.Position = [0.05 0.18 0.93 0.22];

print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_Net_Runtime_All.png')));
plot_num = plot_num + 1;

%% genreate a table of average run times for each illumination/exposure level
fprintf('\n\n');

caption = 'Average DFD-Net Run Time for the Middlebury College Stereo Vision Dataset';
write_latex_table_head(caption, 'tbl:dnn_runtime', 'l|c|c|c|')
fprintf('\\cline{2-4}\n');
fprintf('{} & \\multicolumn{3}{|c|}{\\textbf{Run Time (s)}} \\\\ \\cline{1-4}\n');
fprintf('\\multicolumn{1}{|c|}{\\textbf{Lighting}} & \\textbf{Exposure 0} & \\textbf{Exposure 1} & \\textbf{Exposure 2} \\\\ \\hline \n');
for idx=1:3
    fprintf('\\multicolumn{1}{|l|}{\\textbf{Illumination %d}} & %2.3f & %2.3f & %2.3f \\\\ \\hline \n', idx, avg_runtime(idx,1), avg_runtime(idx,2), avg_runtime(idx,3));
end
write_latex_table_end;

fprintf('\n\n');

fprintf('Average Runtime (s): %2.4f\n\n',mean(data(:,irun)));

%% Sort the data in terms of best to worst from SSIM, NRMSE, NMAE

[data_sort, s_i] = sortrows(data,[iNRMSE,iNMAE],'ascend');

name_sort = {data_name{s_i}}';

img_view = {'Left','Right'};

% print out the table in a data format to copy into latex :-(
caption = 'Top 5 and Bottom 5 DFD-Net Performance Results for the Middlebury College Stereo Vision Dataset';
write_latex_table_head(caption, 'tbl:dnn_top_perf1', '|l|c|c|c|c|c|c|')
fprintf('\\textbf{Name} & \\textbf{View} & \\textbf{Illumination} & \\textbf{Exposure} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\thickhline{2pt} \n');
fprintf('\\multicolumn{7}{|c|}{\\textbf{Top 5}} \\\\ \\cline{1-7}\n');
for idx=1:5
    fprintf('%s & %s & %d & %d & %2.5f & %2.5f & %2.5f \\\\ \\cline{1-7} \n', name_sort{idx}, img_view{data_sort(idx,3)}, data_sort(idx,1), data_sort(idx,2)-1, data_sort(idx,[5,4,6]));
end

fprintf('\\multicolumn{7}{|c|}{\\textbf{Bottom 5}} \\\\ \\cline{1-7}\n');
for idx=length(data)-4:length(data)
    fprintf('%s & %s & %d & %d & %2.5f & %2.5f & %2.5f \\\\ \\cline{1-7} \n', name_sort{idx}, img_view{data_sort(idx,3)}, data_sort(idx,1), data_sort(idx,2)-1, data_sort(idx,[5,4,6]));
end

write_latex_table_end;

fprintf('\n\n');

% write the full data to a table
caption = 'DFD-Net Performance Results for the Middlebury College Stereo Vision Dataset';
write_latex_table_head(caption, 'tbl:dnn_top_perf2', '|l|c|c|c|c|c|c|')
fprintf('\\textbf{Name} & \\textbf{View} & \\textbf{Illumination} & \\textbf{Exposure} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\hline \n');
for idx=1:length(data_sort)
    fprintf('%s & %s & %d & %d & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', name_sort{idx}, img_view{data_sort(idx,3)}, data_sort(idx,1), data_sort(idx,2)-1, data_sort(idx,[5,4,6]));
end

write_latex_table_end;

fprintf('\n\n');

%% print out a table that shows the mean and the variance for each illumination level

for edx=1:3
    NRMSE_mean(edx,1) = mean(NRMSE(1,:,edx),'omitnan');
    NRMSE_mean(edx,2) = mean(NRMSE(2,:,edx),'omitnan');
    NRMSE_mean(edx,3) = mean(NRMSE(3,:,edx),'omitnan');
    NRMSE_std(edx,1) = std(NRMSE(1,:,edx),'omitnan');
    NRMSE_std(edx,2) = std(NRMSE(2,:,edx),'omitnan');
    NRMSE_std(edx,3) = std(NRMSE(3,:,edx),'omitnan');
    
    NMAE_mean(edx,1) = mean(NMAE(1,:,edx),'omitnan');
    NMAE_mean(edx,2) = mean(NMAE(2,:,edx),'omitnan');
    NMAE_mean(edx,3) = mean(NMAE(3,:,edx),'omitnan');
    NMAE_std(edx,1) = std(NMAE(1,:,edx),'omitnan');
    NMAE_std(edx,2) = std(NMAE(2,:,edx),'omitnan');
    NMAE_std(edx,3) = std(NMAE(3,:,edx),'omitnan');
    
    SSIM_mean(edx,1) = mean(SSIM(1,:,edx),'omitnan');
    SSIM_mean(edx,2) = mean(SSIM(2,:,edx),'omitnan');
    SSIM_mean(edx,3) = mean(SSIM(3,:,edx),'omitnan');
    SSIM_std(edx,1) = std(SSIM(1,:,edx),'omitnan');
    SSIM_std(edx,2) = std(SSIM(2,:,edx),'omitnan');
    SSIM_std(edx,3) = std(SSIM(3,:,edx),'omitnan');
    
    NRMSE_mean_exp(edx) = mean([NRMSE(1,:,edx)';NRMSE(2,:,edx)';NRMSE(3,:,edx)'],'omitnan');
    NMAE_mean_exp(edx) = mean([NMAE(1,:,edx)';NMAE(2,:,edx)';NMAE(3,:,edx)'],'omitnan');
    SSIM_mean_exp(edx) = mean([SSIM(1,:,edx)';SSIM(2,:,edx)';SSIM(3,:,edx)'],'omitnan');
    
    NRMSE_std_exp(edx) = std([NRMSE(1,:,edx)';NRMSE(2,:,edx)';NRMSE(3,:,edx)'],'omitnan');
    NMAE_std_exp(edx) = std([NMAE(1,:,edx)';NMAE(2,:,edx)';NMAE(3,:,edx)'],'omitnan');
    SSIM_std_exp(edx) = std([SSIM(1,:,edx)';SSIM(2,:,edx)';SSIM(3,:,edx)'],'omitnan');
end

fprintf('\n\n');

caption = 'DFD-Net Performance Mean \& Standard Deviation';
write_latex_table_head(caption, 'tbl:dnn_mean1', 'l|l|c|c|c|c|c|c|')
fprintf('\\cline{3-8} \n');
fprintf('\\multicolumn{2}{c|}{} & \\multicolumn{2}{c|}{\\textbf{NRMSE}} & \\multicolumn{2}{c|}{\\textbf{NMAE}} & \\multicolumn{2}{c|}{\\textbf{SSIM}} \\\\ \\cline{3-8} \n');
fprintf('\\multicolumn{2}{c|}{} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Mean} & \\textbf{Std} \\\\ \\hline \n');

% print out the values
for edx=1:3
    for illum=1:3
        fprintf('\\multicolumn{1}{|l}{Exp %d} &  \\multicolumn{1}{|l|}{Illum %d} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\cline{2-8} \n',edx-1,illum,NRMSE_mean(edx,illum),NRMSE_std(edx,illum),NMAE_mean(edx,illum),NMAE_std(edx,illum),SSIM_mean(edx,illum),SSIM_std(edx,illum));   
    end 
    fprintf('\\multicolumn{1}{|l}{} &  \\multicolumn{1}{|l|}{Overall} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\thickhline{2pt} \n',NRMSE_mean_exp(edx),NRMSE_std_exp(edx),NMAE_mean_exp(edx),NMAE_std_exp(edx),SSIM_mean_exp(edx),SSIM_std_exp(edx));
end
write_latex_table_end;

fprintf('\n\n');

% NRMSE_mean_comb = mean(data(:,iNRMSE),'omitnan');
% NRMSE_std_comb = std(data(:,iNRMSE),'omitnan');
% 
% NMAE_mean_comb = mean(data(:,iNMAE),'omitnan');
% NMAE_std_comb = std(data(:,iNMAE),'omitnan');
% 
% SSIM_mean_comb = mean(data(:,iSSIM),'omitnan');
% SSIM_std_comb = std(data(:,iSSIM),'omitnan');
% 
% fprintf('\n\n');
% 
% caption = 'DFD-Net Performance Mean \& Standard Deviation';
% write_latex_table_head(caption, 'tbl:dnn_mean1', 'l|l|c|c|c|c|c|c|')
% fprintf('\\cline{3-8} \n');
% fprintf('\\multicolumn{2}{c|}{} & \\multicolumn{2}{c|}{\\textbf{NRMSE}} & \\multicolumn{2}{c|}{\\textbf{NMAE}} & \\multicolumn{2}{c|}{\\textbf{SSIM}} \\\\ \\cline{3-8} \n');
% fprintf('\\multicolumn{2}{c|}{} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Mean} & \\textbf{Std} \\\\ \\hline \n');
% 
% % print out the values
% for illum=1:3
% 
%     for edx=1:3
%         fprintf('\\multicolumn{1}{|l}{Illum %d} &  \\multicolumn{1}{|l|}{Exp %d} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\cline{2-8} \n',illum,edx-1,NRMSE_mean(illum,edx),NRMSE_std(illum,edx),NMAE_mean(illum,edx),NMAE_std(illum,edx),SSIM_mean(illum,edx),SSIM_std(illum,edx));
%     end
% end
% 
% fprintf('\\multicolumn{2}{|l|}{Overall} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\hline \n',NRMSE_mean_comb,NRMSE_std_comb,NMAE_mean_comb,NMAE_std_comb,SSIM_mean_comb,SSIM_std_comb);
% write_latex_table_end;
% 
% fprintf('\n\n');


%% save the results in a mat file for later processing

disp('Saving variables...');
save_name = fullfile(log_path,strcat('dfd_dnn_',scenario,'_results.mat'));

save(save_name, 'data', 'data_name', 'data_sort', 'name_sort', 's_i', 'dm_file', 'focus_file', 'defocus_file', 'gt_file'); 

%% this section will then start to stitch together the images

data_dir = 'D:\IUPUI\Test_Data\Middlebury_Images_Third';

if(generate_images)
    
    img_space_h = 6;     % number of pixels to put between the images
    img_space_v = 6;     % number of pixels to put between the images

    tb_top = s_i([1:5,length(data)-4:length(data)]);

    % image order: infocus, defocused, groundtruth, result
    combined_img = {};

    img_h = 110;
    img_w = 125;

    for idx=1:length(tb_top)

        i_img = imread(fullfile(data_dir, focus_file{tb_top(idx),1}));
        d_img = imread(fullfile(data_dir, defocus_file{tb_top(idx),1}));
        gt_img = imread(fullfile(data_dir, gt_file{tb_top(idx),1}));
        gt_img = cat(3, gt_img, gt_img, gt_img);
        %dm_file = strsplit(focus_file{tb(idx),1},'\');
        dm_img = imread(fullfile(log_path,dm_file{tb_top(idx),1}));
        if(size(dm_img,3) == 1)
            dm_img = cat(3, dm_img, dm_img, dm_img);
        end

        % resize the images
        i_img = imresize(i_img,[img_h,img_w]);
        d_img = imresize(d_img,[img_h,img_w]);
        gt_img = imresize(gt_img,[img_h,img_w]);
        dm_img = imresize(dm_img,[img_h,img_w]);

        im_size = size(dm_img);
        pad_img = 255*ones([im_size(1) img_space_h im_size(3)]);

        i_img = i_img(1:im_size(1), 1:im_size(2), :);
        d_img = d_img(1:im_size(1), 1:im_size(2), :);
        gt_img = gt_img(1:im_size(1), 1:im_size(2), :);

        combined_img{idx} = [i_img pad_img d_img pad_img gt_img pad_img dm_img];

    %     figure(plot_num)
    %     set(gcf,'position',([100,100,1000,300]),'color','w')
    %     image(combined_img{idx});
    %     axis off

        save_file = strcat(data_name{tb_top(idx)}, '_I', num2str(data(tb_top(idx),1)), '_E', num2str(data(tb_top(idx),2)-1),'_V',num2str(data(tb_top(idx),3)),'_',num2str(idx),'_combined.png');   
        imwrite(combined_img{idx},fullfile(log_path,save_file));

    end

    % merge the images together to create the combined top/bottom 5 image
    cmb_w = size(combined_img{1},2);

    pad1 = 255*ones([img_space_v 2*(cmb_w+img_space_h) 3]);
    pad2 = 255*ones([img_h 2*img_space_h 3]);
    full_img = [];

    for idx=1:4

        full_img = [full_img; combined_img{idx} pad2 combined_img{idx+5}];
        full_img = [full_img; pad1];
    end
    full_img = [full_img; combined_img{5} pad2 combined_img{10}];

    save_file = strcat('dfd_dnn_comb_tb_mb.png');
    imwrite(full_img,fullfile(log_path,save_file));

end


%% this section will generate images of all the results in an Nx4 pattern

base_name = 'dfd_mb_pso_all_';

img_row = 7;
clr_select = 1;

img_h = 90;
img_w = 100;

img_space_h = 5;     % number of pixels to put between the images
img_space_v = 5;     % number of pixels to put between the images

num_images = size(s_i, 1);

num_img_grp = floor(num_images/img_row);

num_img_rem = num_images - num_img_grp*img_row;

if(generate_images)
    
    index = [1:1:img_row];
    
    for jdx=1:num_img_grp
    
        img_grp = s_i(index);

        % image order: infocus, defocused, groundtruth, result
        combined_img = {};   
               
        for idx=1:length(img_grp)

            i_img = imread(fullfile(data_dir, focus_file{img_grp(idx),1}));
            d_img = imread(fullfile(data_dir, defocus_file{img_grp(idx),1}));

            gt_img = imread(fullfile(data_dir, gt_file{img_grp(idx),1}));

            dm_img = imread(fullfile(log_path,dm_file{img_grp(idx),1}));
           
            if(clr_select == 2)
                clr_name = 'jet';
                dm_img = 255*ind2rgb(dm_img, jet(256));               
                gt_img = 255*ind2rgb(gt_img, jet(256));
            elseif(clr_select == 1)
                clr_name = 'gry';
                if(size(dm_img,3) == 1)
                    dm_img = cat(3, dm_img, dm_img, dm_img);
                end
                
                if(size(gt_img,3) == 1)
                    gt_img = cat(3, gt_img, gt_img, gt_img);
                end                

            else
                clr_name = '';
            end

            % resize the images
            i_img = imresize(i_img,[img_h,img_w]);
            d_img = imresize(d_img,[img_h,img_w]);
            gt_img = imresize(gt_img,[img_h,img_w]);
            dm_img = imresize(dm_img,[img_h,img_w]);

            im_size = size(dm_img);
            pad_img = 255*ones([im_size(1) img_space_h im_size(3)]);

%             i_img = i_img(1:im_size(1), 1:im_size(2), :);
%             d_img = d_img(1:im_size(1), 1:im_size(2), :);
%             gt_img = gt_img(1:im_size(1), 1:im_size(2), :);

            combined_img{idx} = [i_img pad_img d_img pad_img gt_img pad_img dm_img];

        end

        % merge the images together to create the combined top/bottom 5 image
        cmb_w = size(combined_img{1},2);

        pad1 = 255*ones([img_space_v cmb_w 3]);
        
        full_img = [];

        for idx=1:img_row-1
            full_img = [full_img; combined_img{idx}];
            full_img = [full_img; pad1];
        end
        full_img = [full_img; combined_img{img_row}];
        
        save_file = strcat(base_name, clr_name, '_pt', num2str(jdx,'%02d'), '.jpg');
        imwrite(full_img, fullfile(log_path,save_file));
        
        index = index + img_row;
    
    end
    
    index = index(index<=num_images);
    
    img_grp = s_i(index);

    % image order: infocus, defocused, groundtruth, result
    combined_img = {};   

    for idx=1:length(img_grp)

        i_img = imread(fullfile(data_dir, focus_file{img_grp(idx),1}));
        d_img = imread(fullfile(data_dir, defocus_file{img_grp(idx),1}));

        gt_img = imread(fullfile(data_dir, gt_file{img_grp(idx),1}));

        dm_img = imread(fullfile(log_path,dm_file{img_grp(idx),1}));

        if(clr_select == 2)
            clr_name = 'jet';
            dm_img = 255*ind2rgb(dm_img, jet(256));               
            gt_img = 255*ind2rgb(gt_img, jet(256));
        elseif(clr_select == 1)
            clr_name = 'gry';
            if(size(dm_img,3) == 1)
                dm_img = cat(3, dm_img, dm_img, dm_img);
            end

            if(size(gt_img,3) == 1)
                gt_img = cat(3, gt_img, gt_img, gt_img);
            end                

        else
            clr_name = '';
        end

        % resize the images
        i_img = imresize(i_img,[img_h,img_w]);
        d_img = imresize(d_img,[img_h,img_w]);
        gt_img = imresize(gt_img,[img_h,img_w]);
        dm_img = imresize(dm_img,[img_h,img_w]);

        im_size = size(dm_img);
        pad_img = 255*ones([im_size(1) img_space_h im_size(3)]);

%             i_img = i_img(1:im_size(1), 1:im_size(2), :);
%             d_img = d_img(1:im_size(1), 1:im_size(2), :);
%             gt_img = gt_img(1:im_size(1), 1:im_size(2), :);

        combined_img{idx} = [i_img pad_img d_img pad_img gt_img pad_img dm_img];

    end

    if(~isempty(combined_img))
        % merge the images together to create the combined top/bottom 5 image
        cmb_w = size(combined_img{1},2);

        pad1 = 255*ones([img_space_v cmb_w 3]);

        full_img = [];

        for idx=1:numel(combined_img)-1

            full_img = [full_img; combined_img{idx}];
            full_img = [full_img; pad1];
        end
        full_img = [full_img; combined_img{numel(combined_img)}];

        save_file = strcat(base_name, clr_name, '_pt', num2str(num_img_grp+1, '%02d'), '.jpg');
        imwrite(full_img, fullfile(log_path,save_file));
    end
end