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

cm = [1,0,0; 0,1,0; 0,0,1; 1,0.5,0; 1,1,0; 0.6,0,1; 0,1,1; 1,0,1];

generate_images = false;

if(ispc)
    data_directory = 'D:\IUPUI\Test_Data\';
else
    data_directory = '/home/owner/DfD/data/';
end

%% get the file to scan
startpath = 'D:\IUPUI\PhD\Results\dfd_dnn_rw\';
[log_file, log_path] = uigetfile(file_filter, 'Select DFD-Net Log File', startpath);
if(log_path == 0)
    return;
end

commandwindow;


%% Start reading the file

data_name = {};
focus_file = {};
defocus_file = {};
gt_file = {};
dm_file = {};
data = [];
img_size = [];

% open the file
file_id = fopen(fullfile(log_path,log_file),'r');

index = 1;

while ~feof(file_id)
    
    tline = fgetl(file_id);
    
    % look for the separator
    if(strcmp(tline,sep))
        tline = fgetl(file_id);
        
        % becuase there are many separators in the file llok for the
        % 'Depthmap' word
        if(~isempty(tline) && strcmp(tline(1:8),'Depthmap'))
            
            % run time
            r_line = regexp(tline,': ','split');
            r_line2 = regexp(r_line{2},' ', 'split');
            run_time = str2double(r_line2{1});
            
            % image size
            tline = fgetl(file_id);
            i_line = regexp(tline, ': ', 'split');
            i_line = regexp(i_line{2}, ' x ', 'split');
            img_size(index,:) = [str2double(i_line{1}), str2double(i_line{2})];            
            
            % focus file name
            tline = fgetl(file_id);
            f_line = regexp(tline,'Focus File: ', 'split');
            folders = regexp(f_line{2}, '/', 'split');
            focus_file{index,1} = fullfile(folders{end-4:end});
            
            % get the data name
            data_name{index,1} = folders{numel(folders)-3};
            
            % get the exposure time
            exp_tmp = regexp(folders{end-1}, '_', 'split');
            exp_time = str2double(exp_tmp{2});
            
            % defocus file name
            tline = fgetl(file_id);
            f_line = regexp(tline, 'Defocus File: ', 'split');
            folders = regexp(f_line{2}, '/', 'split');
            defocus_file{index,1} = fullfile(folders{end-4:end});
            
            % ground truth file 
            gt_file{index,1} = fullfile(folders{end-4:end-3},'lidar/lidar_rng_left_00000_8bit_x6.png');
            
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
            
            tline = fgetl(file_id);
            r2 = regexp(tline,': ','split');
            SILOG = str2double(r2{2});            
            
            data(index,:) = [exp_time, NMAE, NRMSE, SSIM, SILOG, run_time];
           
            index = index + 1;            
        end
        
    end
end

fclose(file_id);

clear exp_time NMAE NRMSE SSIM SILOG run_time tline f_line i_line r_line r_line2 r2 folders

fprintf('Number of images to process: %04d\n\n',size(data,1));

%% start processing the data

unique_data_name = unique(data_name);

%indexes for the data
iEXP = 1;
iNMAE = 2;
iNRMSE = 3;
iSSIM = 4;
iSILOG = 5;
irun = 6;
scenario = 'test';

% get the exposure time min, max
min_exp = min(data(:,1));
max_exp = max(data(:,1));
exp_step = 10;

unique_exp_time = num2str([min_exp:exp_step:max_exp]');

% get the mean results for each exposure time
index = 1;
exp_data = [];
exp_nrmse = [];
exp_nmae = [];
exp_ssim = [];
for exp_idx=min_exp:exp_step:max_exp
     
    exp_data(:,:,index) = data(data(:,1) == exp_idx,:);   
    mean_exp_data(index,:) = [exp_idx, mean(exp_data(:,2:end,index),1)];
    
    exp_nrmse(index,:) = exp_data(:,iNRMSE,index);
    exp_nmae(index,:) = exp_data(:,iNMAE,index);
    exp_ssim(index,:) = exp_data(:,iSSIM,index);   
    
    index = index + 1;
    
end

exp_nrmse_std = std(exp_nrmse,0,2);
exp_nmae_std = std(exp_nmae,0,2);
exp_ssim_std = std(exp_ssim,0,2);

%% calc the min, mean and max stats and then print them out

nrmse_stats = [min(data(:,iNRMSE)), mean(data(:,iNRMSE)), max(data(:,iNRMSE)), std(data(:,iNRMSE))];
nmae_stats = [min(data(:,iNMAE)), mean(data(:,iNMAE)), max(data(:,iNMAE)), std(data(:,iNMAE))];
ssim_stats = [min(data(:,iSSIM)), mean(data(:,iSSIM)), max(data(:,iSSIM)), std(data(:,iSSIM))];
stats_name = {'Minimum', 'Mean', 'Maximum', 'Std Deviation'};

fprintf('\n');
caption = 'DfD-Net Real World Dataset Overall Test Performance Results';
write_latex_table_head(caption, 'tbl:dfd_rw_stats1', 'l|c|c|c|')
fprintf('\\cline{2-4}\n');
fprintf('\\multicolumn{1}{c|}{} & {\\textbf{NRMSE}} & {\\textbf{NMAE}} & {\\textbf{SSIM}} \\\\ \\hline \n');
% print out the values
for idx=1:numel(stats_name)    
    fprintf('\\multicolumn{1}{|l|}{%s} & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', stats_name{idx}, nrmse_stats(idx), nmae_stats(idx), ssim_stats(idx));     
end
write_latex_table_end;
fprintf('\n\n');


%% print out the data in latex form
% Average Exposure Time Performance Results

caption = 'DfD-Net Real World Dataset Average Exposure Time Performance Results';
write_latex_table_head(caption, 'tbl:dfd_rw_exp1', '|c|c|c|c|c|c|c|')
%fprintf('{\\textbf{Exposure Time (ms)}} & {\\textbf{NMAE}} & {\\textbf{NRMSE}} & {\\textbf{SSIM}} \\\\ \\hline \n');

fprintf('\\multicolumn{1}{|c}{\\textbf{Exposure}} & \\multicolumn{2}{|c|}{\\textbf{NMAE}} & \\multicolumn{2}{|c|}{\\textbf{NRMSE}} & \\multicolumn{2}{|c|}{\\textbf{SSIM}} \\\\ \\cline{2-7} \n');
fprintf('{\\textbf{Time (ms)}} & {\\textbf{Mean}} & {\\textbf{Std}} & {\\textbf{Mean}} & {\\textbf{Std}} & {\\textbf{Mean}} & {\\textbf{Std}} \\\\ \\hline \n');

% print out the values
for idx=1:size(mean_exp_data, 1)
    
    %fprintf('%d & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', mean_exp_data(idx,1:4));   
    fprintf('%d & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', mean_exp_data(idx,1), mean_exp_data(idx,2), exp_nrmse_std(idx), mean_exp_data(idx,3), exp_nmae_std(idx), mean_exp_data(idx,4), exp_ssim_std(idx));
    
end
write_latex_table_end;

fprintf('\n\n');


%% get the average results per scene

num_scene = numel(unique_data_name);
mean_scenario_data = zeros(num_scene,6);

for idx=1:num_scene
    
    index = find(strcmp(data_name, unique_data_name{idx}));
    
    mean_scenario_data(idx,:) = mean(data(index,:),1);

end


%% plot the average results per scene - All exposure levels combined

figure(plot_num);
set(gcf,'position',([50,50,1400,700]),'color','w')

subplot(3,1,1)
hold on
box on
grid on

b1 = bar([1:num_scene], mean_scenario_data(:,iNRMSE));
b1.BarWidth = 0.6;
b1(1).FaceColor = 'b';

set(gca,'fontweight','bold','FontSize', 13)

% X-Axis
xlim([0,num_scene+1]);
xticks(1:num_scene);
xticklabels([]);

% Y-Axis
plt_step = 33.33333;
nrmse_plt_max = max(mean_scenario_data(:,iNRMSE));
nrmse_plt_max = ceil(nrmse_plt_max*plt_step)/plt_step;
nrmse_plt_min = 0.00;

ylim([nrmse_plt_min nrmse_plt_max]);
yticks(nrmse_plt_min:(1/plt_step):nrmse_plt_max);
ytickformat('%1.2f');
ylabel('NRMSE','fontweight','bold','FontSize', 13)

title(strcat('DfD-Net Synthetically Blurred Real World Performance Results'),'fontweight','bold','FontSize', 16)
ax = gca;
ax.Position = [0.055 0.695 0.93 0.25];

%-------------------------------------------------------------------------
subplot(3,1,2)
hold on
box on
grid on

b1=bar([1:num_scene], mean_scenario_data(:,iNMAE),'b');
b1.BarWidth = 0.6;
b1(1).FaceColor = 'b';

set(gca,'fontweight','bold','FontSize', 13)

% X-Axis
xlim([0,num_scene+1]);
xticks(1:num_scene);
xticklabels([]);

% Y-Axis
plt_step = 50;
nmae_plt_max = max(mean_scenario_data(:,iNMAE));
nmae_plt_max = ceil(nmae_plt_max*plt_step)/plt_step;
nmae_plt_min = 0;

ylim([nmae_plt_min nmae_plt_max]);
yticks(nmae_plt_min:(1/plt_step):nmae_plt_max);
ytickformat('%1.2f');
ylabel('NMAE','fontweight','bold','FontSize', 13)

ax = gca;
ax.Position = [0.055 0.405 0.93 0.25];
    
%-------------------------------------------------------------------------
subplot(3,1,3)
hold on
box on
grid on

b1 = bar([1:num_scene], mean_scenario_data(:,iSSIM),'b');
b1.BarWidth = 0.6;
b1(1).FaceColor = 'b';

set(gca, 'fontweight', 'bold', 'FontSize', 13)

% X-Axis
xlim([0,num_scene+1]);
xticks(1:num_scene);
xticklabels(unique_data_name);
xtickangle(90);
xlabel('Scene', 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
plt_step = 4;
ssim_plt_max = 1;
ssim_plt_min = 0;

ylim([ssim_plt_min ssim_plt_max]);
yticks((ssim_plt_min:(1/plt_step):ssim_plt_max));
ytickformat('%1.2f');
ylabel('SSIM','fontweight','bold','FontSize', 13)

ax = gca;
ax.Position = [0.055 0.115 0.93 0.25];

if(generate_images)
    print(plot_num, '-dpng', fullfile(log_path,strcat('dfd_rw_perf1.png')));
end

plot_num = plot_num + 1;

%% plot the average results for each scene - all exposure times combined and add the mean and std dev

sig_m = 2;

% get the lower limits
limits_test(1,1) = max(nmae_stats(2) - sig_m*nmae_stats(4), 0);
limits_test(2,1) = max(nrmse_stats(2) - sig_m*nrmse_stats(4), 0);
limits_test(3,1) = max(ssim_stats(2) - sig_m*ssim_stats(4), 0);
limits_test(4,1) = max(nmae_stats(2) - 1*nmae_stats(4), 0);
limits_test(5,1) = max(nrmse_stats(2) - 1*nrmse_stats(4), 0);
limits_test(6,1) = max(ssim_stats(2) - 1*ssim_stats(4), 0);

% get the upper limits
limits_test(1,2) = nmae_stats(2) + sig_m*nmae_stats(4);
limits_test(2,2) = nrmse_stats(2) + sig_m*nrmse_stats(4);
limits_test(3,2) = ssim_stats(2) + sig_m*ssim_stats(4);
limits_test(4,2) = nmae_stats(2) + 1*nmae_stats(4);
limits_test(5,2) = nrmse_stats(2) + 1*nrmse_stats(4);
limits_test(6,2) = ssim_stats(2) + 1*ssim_stats(4);


figure(plot_num)
set(gcf,'position',([50,50,1300,700]),'color','w')
%set(gcf,'position',([50,50,800,500]),'color','w')

% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(2,2); num_scene+2,limits_test(2,2); num_scene+2,limits_test(2,1); -1,limits_test(2,1)];
a2 = [-1,limits_test(5,2); num_scene+2,limits_test(5,2); num_scene+2,limits_test(5,1); -1,limits_test(5,1)];

patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
patch('Faces',[1 2 3 4], 'Vertices',a2,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);

% a1=area([-1,num_scene+2],[limits_test(2,2), limits_test(2,2)], 'BaseValue', limits_test(2,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth', 1);
% a2=area([-1,num_scene+2],[limits_test(5,2), limits_test(5,2)], 'BaseValue', limits_test(5,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth', 1);

%plot the mean line - maybe
plot([-1,num_scene+2],[nrmse_stats(2), nrmse_stats(2)], '--b', 'LineWidth', 1);

% plot the min line
plot([-1,num_scene+2],[limits_test(2,1), limits_test(2,1)], '-b', 'LineWidth', 1);
plot([-1,num_scene+2],[limits_test(5,1), limits_test(5,1)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_scene], mean_scenario_data(:,iNRMSE),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_scene+1]);
xticks([0:num_scene+1]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xticklabels([]);

% Y-Axis
%nrmse_plt_max = ceil(limits_test(2,2)*80)/80;
%nrmse_plt_min = floor(limits_test(2,1)*80)/80;
plt_step = 33.33333;
nrmse_plt_max = max(max(mean_scenario_data(:,iNRMSE)),limits_test(2,2));
nrmse_plt_max = ceil(nrmse_plt_max*plt_step)/plt_step;
nrmse_plt_min = 0.00;

ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:(1/plt_step):nrmse_plt_max]);
ytickformat('%1.2f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

title('DfD-Net Real World Scene Average Test Results','fontweight','bold','FontSize', 16)

ax = gca;
% ax.XAxis.MinorTickValues = 1:1:num_scene;
ax.Position = [0.065 0.700 0.92 0.24];
%ax.Position = [0.10 0.710 0.88 0.23];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(1,2); num_scene+2,limits_test(1,2); num_scene+2,limits_test(1,1); -1,limits_test(1,1)];
a2 = [-1,limits_test(4,2); num_scene+2,limits_test(4,2); num_scene+2,limits_test(4,1); -1,limits_test(4,1)];

patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
patch('Faces',[1 2 3 4], 'Vertices',a2,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%area([-1,num_scene+2],[limits_test(1,2), limits_test(1,2)], 'BaseValue', limits_test(1,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_scene+2],[nmae_stats(2), nmae_stats(2)], '--b', 'LineWidth', 1);
plot([-1,num_scene+2],[limits_test(1,1), limits_test(1,1)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_scene], mean_scenario_data(:,iNMAE), 15, 'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_scene+1]);
xticks([0:num_scene+1]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xticklabels([]);

% Y-Axis
%nmae_plt_max = ceil(limits_test(1,2)*200)/200;
%nmae_plt_min = floor(limits_test(1,1)*200)/200;
plt_step = 50;
nmae_plt_min = 0.00;
nmae_plt_max = max(max(mean_scenario_data(:,iNMAE)),limits_test(1,2));
nmae_plt_max = ceil(nmae_plt_max*plt_step)/plt_step;

ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:(1/plt_step):nmae_plt_max]);
ytickformat('%1.2f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
% ax.XAxis.MinorTickValues = 1:1:num_scene;
ax.Position = [0.065 0.425 0.92 0.24];
%ax.Position = [0.10 0.44 0.88 0.23];

%--------------------------------------
subplot(3,1,3);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(3,2); num_scene+2,limits_test(3,2); num_scene+2,limits_test(3,1); -1,limits_test(3,1)];
a2 = [-1,limits_test(6,2); num_scene+2,limits_test(6,2); num_scene+2,limits_test(6,1); -1,limits_test(6,1)];

p1=patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
patch('Faces',[1 2 3 4], 'Vertices',a2,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%p1 = area([-1,num_scene+2],[limits_test(3,2), limits_test(3,2)], 'BaseValue', limits_test(3,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1);

%plot the mean line - maybe
p2 = plot([-1,num_scene+2],[ssim_stats(2), ssim_stats(2)], '--b', 'LineWidth', 1);
plot([-1,num_scene+2],[limits_test(3,1), limits_test(3,1)], '-b', 'LineWidth', 1);

% plot the test points
p3 = scatter([1:num_scene], mean_scenario_data(:,iSSIM),15,'filled','k');
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
% xlabel_name = {''};
% for idx=1:num_trials
%     xlabel_name{idx+1,1} = num2str(idx);
% end
xlim([0 num_scene+1]);
xticks([1:num_scene]);
xticklabels(unique_data_name);
xtickangle(90);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');

xlabel(strcat('Scene'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
%ssim_plt_max = ceil(limits_test(3,2)*80)/80;
%ssim_plt_min = floor(limits_test(3,1)*80)/80;

plt_step = 10;
ssim_plt_min = min(min(mean_scenario_data(:,iSSIM)),limits_test(3,1));
ssim_plt_min = floor(ssim_plt_min*plt_step)/plt_step;
ssim_plt_max = max(max(mean_scenario_data(:,iSSIM)),limits_test(3,2));
ssim_plt_max = ceil(ssim_plt_max*plt_step)/plt_step;

% ssim_plt_max = 1.1;
% ssim_plt_min = 0.85;
ylim([ssim_plt_min ssim_plt_max]);
yticks([ssim_plt_min:(1/plt_step):ssim_plt_max]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
ytickformat('%1.2f')
ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

%legend([p1, p2, p3], strcat(32,'{\pm}',num2str(sig_m),'{\sigma} Test Bounds'), 'Metric Mean', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')
legend([p1, p2, p3], {strcat(32,'Test Bounds',32), 'Metric Mean', 'Test Results'}, 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
% ax.XAxis.MinorTickValues = 1:1:num_trials;
ax.Position = [0.065 0.145 0.92 0.24];
%ax.Position = [0.1 0.17 0.88 0.23];

if(generate_images)
    print(plot_num, '-dpng', fullfile(log_path,strcat('dfd_rw_perf2.png')));
end

plot_num = plot_num + 1;  

%% Sort the data in terms of best to worst from SSIM, NRMSE, NMAE
% Top 5 and Bottom 5 DfD-Net Performance Results for the Real World Dataset

[data_sort, s_i] = sortrows(data,[iNRMSE,iNMAE],'ascend');

name_sort = {data_name{s_i}}';

% print out the table in a data format to copy into latex :-(
caption = 'Top 5 and Bottom 5 DfD-Net Performance Results for the Real World Dataset';
write_latex_table_head(caption, 'tbl:dnn_rw4_perf1', '|l|c|c|c|c|')
fprintf('\\textbf{Name} & \\begin{tabular}{c}\\textbf{Exposure} \\\\[-1mm] \\textbf{Time (ms)}\\end{tabular} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\hline \n');
fprintf('\\multicolumn{5}{|c|}{\\textbf{Top 5}} \\\\ \\cline{1-5}\n');
for idx=1:5
    fprintf('%s & %d & %2.5f & %2.5f & %2.5f \\\\ \\cline{1-5} \n', name_sort{idx}, data_sort(idx,iEXP), data_sort(idx,[iNRMSE, iNMAE, iSSIM]));
end

fprintf('\\multicolumn{5}{|c|}{\\textbf{Bottom 5}} \\\\ \\cline{1-5}\n');
for idx=length(data)-4:length(data)
    fprintf('%s & %d & %2.4f & %2.4f & %2.4f \\\\ \\cline{1-5} \n', name_sort{idx}, data_sort(idx,iEXP), data_sort(idx,[iNRMSE, iNMAE, iSSIM]));
end

write_latex_table_end;

%% Write all of the results out to the screen 
% All DfD-Net Performance Results for the Real World Dataset

fprintf('\n\n');
% print out the table in a data format to copy into latex :-(
caption = 'DfD-Net Performance Results for the Real World Dataset';
%write_latex_table_head(caption, 'tbl:dnn_rw4_perf_all1', '|l|c|c|c|c|')
%fprintf('\\textbf{Name} & \\begin{tabular}{c}\\textbf{Exposure} \\\\[-1mm] \\textbf{Time (ms)}\\end{tabular} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\hline \n');
%fprintf('\\multicolumn{5}{|c|}{\\textbf{Top 5}} \\\\ \\cline{1-5}\n');
write_latex_table_head(caption, 'tbl:dnn_rw4_perf_all1', '|c|c|c|c|c|c|')
fprintf('\\textbf{Number} & \\textbf{Name} & \\begin{tabular}{c}\\textbf{Exposure} \\\\[-1mm] \\textbf{Time (ms)}\\end{tabular} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\hline \n');
fprintf('\\multicolumn{6}{|c|}{\\textbf{Top 5}} \\\\ \\cline{1-6}\n');
for idx=1:length(data_sort)
    %fprintf('%s & %d & %2.4f & %2.4f & %2.4f \\\\ \\cline{1-5} \n', name_sort{idx}, data_sort(idx,iEXP), data_sort(idx,[iNRMSE, iNMAE, iSSIM]));
    fprintf('%04d, %s & %d & %2.4f & %2.4f & %2.4f \\\\ \\cline{1-6} \n', idx, name_sort{idx}, data_sort(idx,iEXP), data_sort(idx,[iNRMSE, iNMAE, iSSIM]));
end
write_latex_table_end;
fprintf('\n\n');


%% Plot the metrics on a bar3 chart to show the scene vs exposure time results

cm_max = 100;
cm = hot(cm_max);

figure(plot_num);
set(gcf,'position',([20,50,1500,700]),'color','w')

subplot(3,1,1)
hold on;
box on;
grid on;

b = bar3([min_exp:exp_step:max_exp], exp_nrmse);

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

%plt_step = 33.333333;
plt_step = 50;
nrmse_plt_max = ceil(max(exp_nrmse(:))*plt_step)/plt_step;

colormap(cm);
cb = colorbar('fontweight','bold','FontSize', 13, 'Location', 'eastoutside');
[cb.Ticks , cb.TickLabels] = calc_limits(0, nrmse_plt_max, plt_step, '%1.2f');
% cb.Ticks = [0:0.03:ceil(nrmse_max*20)/20];
% cb.TickLabels = num2str(cb.Ticks','%1.2f');
cb.Limits = [0 nrmse_plt_max];

set(gca, 'fontweight','bold','FontSize', 12);

% X-Axis
xlim([0,num_scene+1]);
xticks(1:num_scene);
xticklabels([]);
% xticklabels(unique_data_name);
% xtickangle(-40);
% xlabel('Scene', 'fontweight', 'bold', 'FontSize', 13);

% Y
% unique_exp_time
ylim([(min_exp-5) (max_exp+5)]);
yticks(min_exp:exp_step:max_exp);
ylabel('Exposure (ms)','fontweight','bold','FontSize', 13)

% Z
zlim([0 nrmse_plt_max]);
zticks([]);
zlabel('NRMSE','fontweight','bold','FontSize', 13)
zticklabels([]);

title('Exposure Time Performance Comparison','fontweight','bold','FontSize', 16)

view(170, 45);
ax = gca;
ax.XDir = 'reverse';

ax.YLabel.Rotation = -25;

%ax.Position = [0.05 0.23 0.875 0.7];
ax.Position = [0.055 0.72 0.87 0.23];

subplot(3,1,2)
hold on;
box on;
grid on;

b = bar3([min_exp:exp_step:max_exp], exp_nmae);

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

plt_step = 50;
nmae_plt_max = ceil(max(exp_nmae(:))*plt_step)/plt_step;

colormap(cm);
cb = colorbar('fontweight','bold','FontSize', 13, 'Location', 'eastoutside');
[cb.Ticks , cb.TickLabels] = calc_limits(0, nmae_plt_max, plt_step, '%1.2f');
% cb.Ticks = [0:0.02:ceil(nmae_max*10)/10];
% cb.TickLabels = num2str(cb.Ticks','%1.2f');
cb.Limits = [0 nmae_plt_max];

set(gca, 'fontweight','bold','FontSize', 12);

% X-Axis
xlim([0,num_scene+1]);
xticks(1:num_scene);
xticklabels([]);
% xticklabels(unique_data_name);
% xtickangle(-40);
% xlabel('Scene', 'fontweight', 'bold', 'FontSize', 13);

% Y
% unique_exp_time
ylim([(min_exp-5) (max_exp+5)]);
yticks(min_exp:exp_step:max_exp);
ylabel('Exposure (ms)','fontweight','bold','FontSize', 13)

% Z
zlim([0 nmae_plt_max]);
zticks([]);
zlabel('NMAE','fontweight','bold','FontSize', 13)
zticklabels([]);

view(170, 45);
ax = gca;
ax.XDir = 'reverse';

ax.YLabel.Rotation = -25;

%ax.Position = [0.05 0.23 0.875 0.7];
ax.Position = [0.055 0.42 0.87 0.23];

subplot(3,1,3)
hold on;
box on;
grid on;

b = bar3([min_exp:exp_step:max_exp], 1-exp_ssim);

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

plt_step = 12.5;
ssim_plt_max = ceil(max(1-exp_ssim(:))*plt_step)/plt_step;

colormap(cm);
cb = colorbar('fontweight','bold','FontSize', 13, 'Location', 'eastoutside');
[cb.Ticks , cb.TickLabels] = calc_limits(0, ssim_plt_max, plt_step, '%1.2f');
%cb.Ticks = [0:0.25:1];
%cb.TickLabels = num2str(cb.Ticks','%1.2f');
cb.Limits = [0 ssim_plt_max];

set(gca, 'fontweight','bold','FontSize', 12);

% X-Axis
xlim([0,num_scene+1]);
xticks(1:num_scene);
xticklabels(unique_data_name);
xtickangle(-35);
xlabel('Scene', 'fontweight', 'bold', 'FontSize', 13);

% Y
% unique_exp_time
ylim([(min_exp-5) (max_exp+5)]);
yticks(min_exp:exp_step:max_exp);
ylabel('Exposure (ms)','fontweight','bold','FontSize', 13)

% Z
zlim([0 ssim_plt_max]);
zticks([]);
zlabel('1-SSIM','fontweight','bold','FontSize', 13)
zticklabels([]);

view(170, 45);
ax = gca;
ax.XDir = 'reverse';

ax.YLabel.Rotation = -25;

%ax.Position = [0.05 0.23 0.875 0.7];
ax.Position = [0.055 0.12 0.87 0.23];

if(generate_images)
    print(plot_num, '-dpng', fullfile(log_path,strcat('dfd_rw4_exp_perf.png')));
end 
plot_num = plot_num + 1;

%% save the results in a mat file for later processing

disp('Saving variables...');
save_name = fullfile(log_path,strcat('dfd_dnn_',scenario,'_results.mat'));

save(save_name, 'data', 'data_name', 'data_sort', 'name_sort', 's_i', 'dm_file', 'focus_file', 'defocus_file', 'gt_file'); 


%% this section will then start to stitch together the images

if(generate_images)
    
    img_space_h = 6;     % number of pixels to put between the images
    img_space_v = 6;     % number of pixels to put between the images

    tb_top = s_i([1:5,length(data)-4:length(data)]);

    % image order: infocus, defocused, groundtruth, result
    combined_img = {};

    img_h = 110;
    img_w = 125;
    
    figure(plot_num);
    set(gcf,'position',([100,100,1000,300]),'color','w');
    
    for idx=1:length(tb_top)

        i_img = imread(fullfile(data_directory, focus_file{tb_top(idx),1}));
        d_img = imread(fullfile(data_directory, defocus_file{tb_top(idx),1}));
        gt_img = imread(fullfile(data_directory, gt_file{tb_top(idx),1}));
        %gt_img = cat(3, gt_img, gt_img, gt_img);
        gt_img = 255*ind2rgb(gt_img, jet(256));
        
        dm_img = imread(fullfile(log_path,dm_file{tb_top(idx),1}));
        dm_img = 255*ind2rgb(dm_img, jet(256));
        
%         if(size(dm_img,3) == 1)
%             dm_img = cat(3, dm_img, dm_img, dm_img);
%         end

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
        image(combined_img{idx});
        axis off

        save_file = strcat(data_name{tb_top(idx)}, '_E', num2str(data(tb_top(idx),iEXP)), '_', num2str(idx), '_combined.png');   
        imwrite(combined_img{idx},fullfile(log_path,save_file));

    end
    plot_num = plot_num + 1;

    
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

    figure(plot_num);
    set(gcf,'position',([100,100,1200,600]),'color','w');
    image(full_img);
    axis off
    plot_num = plot_num + 1;
    
    save_file = strcat('dfd_rw4_comb_tb.png');
    imwrite(full_img, fullfile(log_path,save_file));

end

%% this section will generate images of all the results in an Nx4 pattern
invert_gt = true;
base_name = 'dfd_rw4_sb_k01_k03_';

img_row = 5;
clr_select = 2;

img_h = 105;
img_w = 105;

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

            i_img = imread(fullfile(data_directory, focus_file{img_grp(idx),1}));
            d_img = imread(fullfile(data_directory, defocus_file{img_grp(idx),1}));

            gt_img = imread(fullfile(data_directory, gt_file{img_grp(idx),1}));

            dm_img = imread(fullfile(log_path,dm_file{img_grp(idx),1}));
            
            if(invert_gt)
                gt_img = 255 - gt_img;
            end
           
            if(clr_select == 1)
                clr_name = 'jet';
                dm_img = 255*ind2rgb(dm_img, jet(256));               
                gt_img = 255*ind2rgb(gt_img, jet(256));
            elseif(clr_select == 2)
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

%             image(combined_img{idx});
%             axis off
% 
%             save_file = strcat(data_name{tb_top(idx)}, '_E', num2str(data(tb_top(idx),iEXP)), '_', num2str(idx), '_combined.png');   
%             imwrite(combined_img{idx},fullfile(log_path,save_file));

        end

        % merge the images together to create the combined top/bottom 5 image
        cmb_w = size(combined_img{1},2);

        pad1 = 255*ones([img_space_v cmb_w 3]);
        %pad2 = 255*ones([img_h 2*img_space_h 3]);
        
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

        i_img = imread(fullfile(data_directory, focus_file{img_grp(idx),1}));
        d_img = imread(fullfile(data_directory, defocus_file{img_grp(idx),1}));

        gt_img = imread(fullfile(data_directory, gt_file{img_grp(idx),1}));
       
        if(invert_gt)
            gt_img = 255 - gt_img;
        end
        
        dm_img = imread(fullfile(log_path,dm_file{img_grp(idx),1}));

        if(clr_select == 1)
            clr_name = 'jet';
            dm_img = 255*ind2rgb(dm_img, jet(256));               
            gt_img = 255*ind2rgb(gt_img, jet(256));
        elseif(clr_select == 2)
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

