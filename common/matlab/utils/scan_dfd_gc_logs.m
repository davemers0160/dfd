format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ~] = fileparts(full_path);
plot_num = 1;
sep = '------------------------------------------------------------------';

%% get the directory to scan
startpath = 'D:\IUPUI\PhD\Results\dfd_gc';

log_path = uigetdir(startpath,'Select Graph Cuts Result Folder');

if(log_path == 0)
    return;
end

listing = dir(log_path);
listing = listing(3:end);

fprintf('Found %03d listings\n', numel(listing));

%% create a structure to hold the data
data_name = {};
infocus_img_file = {};
defocus_img_file = {};
gt_img_file = {};
dm_img_file = {};
image_files = {};

data = [];
index = 1;

for idx=1:length(listing)
    
    if(listing(idx).isdir)
        
        logfile = dir(strcat(fullfile(listing(idx).folder,listing(idx).name),'\*.txt'));
        
%         [d, d_n, i_f] = read_dfd_gc_log_single(logfile, listing(idx), sep);
%         data(end+1,:) = d;
%         data_name{end+1,1} = d_n;
%         image_files{end+1,1} = i_f;        
        
        [d, d_n, i_f] = read_dfd_gc_log_multi(logfile, listing(idx), sep);
        
%         [data(end+1,:), data_name(end+1,:), image_files(end+1)] = read_dfd_gc_log_multi(logfile, listing(idx), sep);
        data(end+1,:) = d;
        data_name(end+1,:) = d_n;
        image_files(end+1) = i_f;
        
%         if(isempty(logfile))
%            fprintf('No log file found for index: %d\n',idx-1);
%            continue;
%         end
%         
%         dm_file = dir(strcat(fullfile(listing(idx).folder,listing(idx).name),'\*.png'));
%         dm_img_file{index,1} = fullfile(dm_file(2).folder,dm_file(2).name);
%         
%         file_id = fopen(fullfile(logfile(1).folder,logfile(1).name),'r');
%         
%         platform = fgetl(file_id); % read in the fisrt line - platform
%         tline = fgetl(file_id); % read in the second line - focus file name
%         
%         % parse the path into its folder parts
%         folders = regexp(tline,'/','split');
%         
%         infocus_img_file{index,1} = fullfile(folders{numel(folders)-3:numel(folders)});
%         
%         % the structure will be in two parts
%         % the name will be a cell array of string values
%         % the data will be in a numerical array
%         data_name{index,1} = folders{numel(folders)-3};
%         
%         switch folders{numel(folders)-2}
%             
%             case 'Illum1'
%                 I = 1;
%             case 'Illum2'
%                 I = 2;
%             case 'Illum3'
%                 I = 3;
%         end
%         
%         switch folders{numel(folders)-1}
%             case 'Exp0'
%                 E = 1;
%             case 'Exp1'
%                 E = 2;
%             case 'Exp2'
%                 E = 3; 
%         end
%         
%         switch folders{numel(folders)}
%             case 'view1.png'
%                 V = 1;
%             case 'view5.png'
%                 V = 2;
%         end
%         
%         % get the defocused image location
%         tline = fgetl(file_id);
%         folders = regexp(tline,'/','split');
%         defocus_img_file{index,1} = fullfile(folders{numel(folders)-3:numel(folders)});
% 
%         % get the defocused image location
%         tline = fgetl(file_id);
%         folders = regexp(tline,'/','split');
%         gt_img_file{index,1} = fullfile(folders{numel(folders)-1:numel(folders)});
%         
%         run_time = 0;
%         while ~feof(file_id)
% 
%             tline = fgetl(file_id);
%             
%             s_line = regexp(tline,': ','split');
% 
%             if(strcmp(s_line{1}, 'Depth Map Generation Complete (minutes)'))
%                 run_time = str2double(s_line{2});
%             end
%             
%             if(strcmp(tline,sep))
%                 break;
%             end
%         end
% 
%         tline = fgetl(file_id);
%         
%         
%         results = regexp(tline,'\t','split');
%         r2 = regexp(results{1},': ','split');
%         
%         r2_2 = regexp(r2{2},', ', 'split');
%         NMAE = str2double(r2_2{1});
%         NRMSE = str2double(r2_2{2});
%         %r2 = regexp(results{4},': ','split');
%         SSIM = str2double(r2_2{3});
%         
%         %logfile2{index} = logfile(1).name;
% 
%         data(index,:) = [I, E, V, NMAE, NRMSE, SSIM, run_time];
% 
%         fclose(file_id);
%         index = index + 1;
    end
end

%clear NRMSE NMAE SSIM

fprintf('\nParsed %d directories with logfiles.\n\n',size(data,1));


%%
%indexes for the data
iNMAE = 4;
iNRMSE = 5;
iSSIM = 6;
irun = 7;
avg_runtime = [];

% do this to get a baseline set of names for the xticklabel and to fill in
% the blanks for datasets with no data
plot_I = data(data(:,1) == 1,:);
I_name = data_name(data(:,1) == 1,:);
x_name = I_name(plot_I(:,2) == 1,:);

cm = hsv(3);

for illum=1:3
    plot_I = data(data(:,1) == illum,:);
    I_name = data_name(data(:,1) == illum,:);

    plot_E = {};
    E_name = {};
    
    for idx=1:3
        plot_E{idx} = plot_I(plot_I(:,2) == idx,:);
        E_name{idx} = I_name(plot_I(:,2) == idx,:);
    end
 
    avg_runtime(illum,1) = mean(plot_E{1}(:,irun)*60);
    avg_runtime(illum,2) = mean(plot_E{2}(:,irun)*60);
    avg_runtime(illum,3) = mean(plot_E{3}(:,irun)*60);    
      
    for idx=1:numel(x_name)
        if(strcmp(x_name{idx},E_name{1,1}{idx})==0)
            E_name{1,1} = [{E_name{1,1}{1:idx-1}}';{x_name{idx}};{E_name{1,1}{idx:end}}'];
            plot_E{1,1} = [plot_E{1,1}(1:idx-1,:); [illum, 2, nan, nan, nan, nan, nan]; plot_E{1,1}(idx:end,:)];
        end         
        if(strcmp(x_name{idx},E_name{1,2}{idx})==0)
            E_name{1,2} = [{E_name{1,2}{1:idx-1}}';{x_name{idx}};{E_name{1,2}{idx:end}}'];
            plot_E{1,2} = [plot_E{1,2}(1:idx-1,:); [illum, 2, nan, nan, nan, nan, nan]; plot_E{1,2}(idx:end,:)];
        end    
        if(strcmp(x_name{idx},E_name{1,3}{idx})==0)
            E_name{1,3} = [{E_name{1,3}{1:idx-1}}';{x_name{idx}};{E_name{1,3}{idx:end}}'];
            plot_E{1,3} = [plot_E{1,3}(1:idx-1,:); [illum, 3, nan, nan, nan, nan, nan]; plot_E{1,3}(idx:end,:)];
        end
    end
    
    NRMSE(illum,:,1) = plot_E{1}(:,iNRMSE);
    NRMSE(illum,:,2) = plot_E{2}(:,iNRMSE);
    NRMSE(illum,:,3) = plot_E{3}(:,iNRMSE);
    
%     NRMSE_std(illum,1) = std(plot_E{1}(:,iNRMSE),'omitnan');
%     NRMSE_std(illum,2) = std(plot_E{2}(:,iNRMSE),'omitnan');
%     NRMSE_std(illum,3) = std(plot_E{3}(:,iNRMSE),'omitnan');
%     NRMSE_std(illum,4) = std([plot_E{1}(:,iNRMSE);plot_E{2}(:,iNRMSE);plot_E{3}(:,iNRMSE)],'omitnan');

    NMAE(illum,:,1) = plot_E{1}(:,iNMAE);
    NMAE(illum,:,2) = plot_E{2}(:,iNMAE);
    NMAE(illum,:,3) = plot_E{3}(:,iNMAE);
    %NMAE(illum,4) = mean([plot_E{1}(:,iNMAE);plot_E{2}(:,iNMAE);plot_E{3}(:,iNMAE)],'omitnan');
    
%     NMAE_std(illum,1) = std(plot_E{1}(:,iNMAE),'omitnan');
%     NMAE_std(illum,2) = std(plot_E{2}(:,iNMAE),'omitnan');
%     NMAE_std(illum,3) = std(plot_E{3}(:,iNMAE),'omitnan');
%     NMAE_std(illum,4) = std([plot_E{1}(:,iNMAE);plot_E{2}(:,iNMAE);plot_E{3}(:,iNMAE)],'omitnan');
    
    SSIM(illum,:,1) = plot_E{1}(:,iSSIM);
    SSIM(illum,:,2) = plot_E{2}(:,iSSIM);
    SSIM(illum,:,3) = plot_E{3}(:,iSSIM);
    %SSIM(illum,4) = mean([plot_E{1}(:,iSSIM);plot_E{2}(:,iSSIM);plot_E{3}(:,iSSIM)],'omitnan');
    
%     SSIM_std(illum,1) = std(plot_E{1}(:,iSSIM),'omitnan');
%     SSIM_std(illum,2) = std(plot_E{2}(:,iSSIM),'omitnan');
%     SSIM_std(illum,3) = std(plot_E{3}(:,iSSIM),'omitnan');   
%     SSIM_std(illum,4) = std([plot_E{1}(:,iSSIM);plot_E{2}(:,iSSIM);plot_E{3}(:,iSSIM)],'omitnan');
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
%     print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_GC_NRMSE_Results_I',num2str(illum),'.png')));
% 
%     plot_num = plot_num + 1;
% 
%     % plot the NMAE and SSIM
% 
%     nmae_max = max(data(:,iNMAE));
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
%     ylim([0 ceil(nmae_max)])
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
%     print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_GC_NMAE_Results_I',num2str(illum),'.png')));
% 
%     plot_num = plot_num + 1;
%     
    
    run_times{illum,1} = plot_E{1}(:,irun)*60;
    run_times{illum,2} = plot_E{2}(:,irun)*60;
    run_times{illum,3} = plot_E{3}(:,irun)*60;
    
%     % plot the run times for each exposure level
%     figure(plot_num);
%     set(gcf,'position',([100,100,1000,600]),'color','w')
%     hold on
%     box on
%     grid on
%     plot(plot_E{1}(:,irun)*60,'-b')
%     plot(plot_E{2}(:,irun)*60,'--b')
%     plot(plot_E{3}(:,irun)*60,':b')
%     set(gca,'fontweight','bold')
%     xlim([1 length(E_name{1})]);
%     xticks([1:length(E_name{1})]);
%     xticklabels(E_name{1});
%     xtickangle(90);
%     title(strcat('Graph Cuts Runtime - Illumination',32,num2str(illum)))
%     ylabel('Run Time (s)', 'fontweight','bold')
%     legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
%     ax = gca;
%     ax.Position = [0.07 0.23 0.90 0.72];
%  
%     print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_GC_Runtime_I',num2str(illum),'.png')));
%     plot_num = plot_num + 1;
    
    % now plot everything in its own subplot just to see how it looks
    fig = figure(plot_num);
    set(gcf,'position',([50,50,1400,700]),'color','w')
    nrmse_max = max(data(:,iNRMSE));
    
    subplot(3,1,1)
    hold on
    box on
    grid on

    b1=bar([plot_E{1}(:,iNRMSE),plot_E{2}(:,iNRMSE),plot_E{3}(:,iNRMSE)]);
    set(gca,'fontweight','bold','FontSize', 13)
    
%     plot(plot_E{1}(:,iNRMSE),'-b')
%     plot(plot_E{2}(:,iNRMSE),'--b')
%     plot(plot_E{3}(:,iNRMSE),':b')

    % X-Axis
    xlim([0 length(E_name{1})+1]);
    xticks([1:length(E_name{1})]);
    xticklabels([]);
    xtickangle(90);
    
    % Y-Axis
    ylim([0 ceil(10*nrmse_max)/10])
    ylabel('NRMSE','fontweight','bold','FontSize', 13)
    title(strcat('Graph Cuts Performance Results - Illumination',32,num2str(illum)),'fontweight','bold','FontSize', 16)

    b1(1).FaceColor = cm(1,:);
    b1(2).FaceColor = cm(2,:);
    b1(3).FaceColor = cm(3,:);    
    ax = gca;
    ax.Position = [0.05 0.73 0.93 0.22];
    
    nmae_max = max(data(:,iNMAE));    
    subplot(3,1,2)
    hold on
    box on
    grid on

    b2=bar([plot_E{1}(:,iNMAE),plot_E{2}(:,iNMAE),plot_E{3}(:,iNMAE)]);
    set(gca,'fontweight','bold','FontSize', 13)
%     plot(plot_E{1}(:,iNMAE),'-b')
%     plot(plot_E{2}(:,iNMAE),'--b')
%     plot(plot_E{3}(:,iNMAE),':b')

    % X-Axis
    xlim([0 length(E_name{1})+1]);
    xticks([1:length(E_name{1})]);
    xticklabels([]);
    xtickangle(90);
    
    % Y-Axis
    ylim([0 ceil(10*nmae_max)/10])
    ylabel('NMAE','fontweight','bold','FontSize', 13)

    colormap(cm);
    b2(1).FaceColor = cm(1,:);
    b2(2).FaceColor = cm(2,:);
    b2(3).FaceColor = cm(3,:); 
    ax = gca;
    ax.Position = [0.05 0.475 0.93 0.22];

    subplot(3,1,3)
    hold on
    box on
    grid on

    b3=bar([plot_E{1}(:,iSSIM),plot_E{2}(:,iSSIM),plot_E{3}(:,iSSIM)]);
    set(gca,'fontweight','bold','FontSize', 13)
%     plot(plot_E{1}(:,iSSIM),'-b')
%     plot(plot_E{2}(:,iSSIM),'--b')
%     plot(plot_E{3}(:,iSSIM),':b')

    % X-Axis
    xlim([0 length(E_name{1})+1]);
    xticks([1:length(E_name{1})]);
    xticklabels(E_name{1});
    xtickangle(90);
    
    % Y-Axis
    ylim([0 1])
    ylabel('SSIM','fontweight','bold','FontSize', 13)
    yticks((0:0.25:1));
    
    b3(1).FaceColor = cm(1,:);
    b3(2).FaceColor = cm(2,:);
    b3(3).FaceColor = cm(3,:); 
    legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
    ax = gca;
    ax.Position = [0.05 0.22 0.93 0.22];

    print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_GC_Perf_I',num2str(illum),'.png')));

    plot_num = plot_num + 1; 

end

%% generate a runtime graph for each of the illuminations in subplots
illum = 1;

plot_I = data(data(:,1) == illum,:);
I_name = data_name(data(:,1) == illum,:);

E_name = {};
E_name{1} = I_name(plot_I(:,2) == 1,:);

figure(plot_num);
set(gcf,'position',([50,50,1400,700]),'color','w')

subplot(3,1,1)
hold on
box on
grid on
% plot(run_times{illum,1},'-b')
% plot(run_times{illum,2},'--b')
% plot(run_times{illum,3},':b')
b1 = bar([run_times{illum,1},run_times{illum,2},run_times{illum,3}]);
set(gca,'fontweight','bold')
xlim([0 length(E_name{1})+1]);
xticks([1:length(E_name{1})]);
xticklabels([]);
xtickangle(90);
b1(1).FaceColor = cm(1,:);
b1(2).FaceColor = cm(2,:);
b1(3).FaceColor = cm(3,:); 
title(strcat('Graph Cuts Run Time - Illumination',32,num2str(illum)))
ylabel('Run Time (s)', 'fontweight','bold')
%legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.05 0.73 0.93 0.22];

illum = 2;
subplot(3,1,2)
hold on
box on
grid on
% plot(run_times{illum,1},'-b')
% plot(run_times{illum,2},'--b')
% plot(run_times{illum,3},':b')
b2 = bar([run_times{illum,1},run_times{illum,2},run_times{illum,3}]);
set(gca,'fontweight','bold')
xlim([0 length(E_name{1})+1]);
xticks([1:length(E_name{1})]);
xticklabels([]);
xtickangle(90);
b2(1).FaceColor = cm(1,:);
b2(2).FaceColor = cm(2,:);
b2(3).FaceColor = cm(3,:); 
title(strcat('Graph Cuts Run Time - Illumination',32,num2str(illum)))
ylabel('Run Time (s)', 'fontweight','bold')
%legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.05 0.45 0.93 0.22];

illum = 3;
subplot(3,1,3)
hold on
box on
grid on
% plot(run_times{illum,1},'-b')
% plot(run_times{illum,2},'--b')
% plot(run_times{illum,3},':b')
b3 = bar([run_times{illum,1},run_times{illum,2},run_times{illum,3}]);
set(gca,'fontweight','bold')
xlim([0 length(E_name{1})+1]);
xticks([1:length(E_name{1})]);
xticklabels(E_name{1});
xtickangle(90);
b3(1).FaceColor = cm(1,:);
b3(2).FaceColor = cm(2,:);
b3(3).FaceColor = cm(3,:); 
title(strcat('Graph Cuts Run Time - Illumination',32,num2str(illum)))
ylabel('Run Time (s)', 'fontweight','bold')
legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.05 0.18 0.93 0.22];

print(plot_num, '-dpng', fullfile(log_path,strcat('DFD_GC_Runtime_All.png')));
plot_num = plot_num + 1;

%% genreate a table of average run times for each illumination/exposure level
fprintf('\n\n');

caption = 'Average Graph Cut Run Time for the Middlebury College Stereo Vision Dataset';
write_latex_table_head(caption, 'tbl:gc_runtime', 'l|c|c|c|')
fprintf('\\cline{2-4}\n');
fprintf('{} & \\textbf{Exposure 0} & \\textbf{Exposure 1} & \\textbf{Exposure 2} \\\\ \\hline \n');
for idx=1:3
    fprintf('\\multicolumn{1}{|l|}{Illumination %d} & %2.3f & %2.3f & %2.3f \\\\ \\hline \n', idx, avg_runtime(idx,1), avg_runtime(idx,2), avg_runtime(idx,3));
end
write_latex_table_end;

fprintf('\n\n');


%% Sort the data in terms of best to worst from SSIM, NRMSE, NMAE

[data_sort, s_i] = sortrows(data,[iNRMSE,iNMAE,iSSIM],'ascend');

name_sort = {data_name{s_i}}';

img_view = {'Left','Right'};

% print out the table in a data format to copy into latex :-(
caption = 'Top 5 and Bottom 5 Graph Cuts Performance Results for the Middlebury College Stereo Vision Dataset';
write_latex_table_head(caption, 'tbl:gc_perf1', '|l|c|c|c|c|c|c|')
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

% write the full performance results to a table
caption = 'Graph Cuts Performance Results for the Middlebury College Stereo Vision Dataset';
write_latex_table_head(caption, 'tbl:gc_perf2', '|l|c|c|c|c|c|c|')
fprintf('\\textbf{Name} & \\textbf{View} & \\textbf{Illumination} & \\textbf{Exposure} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\hline \n');
for idx=1:length(data)
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

caption = 'Graph Cuts Performance Mean \& Standard Deviation';
write_latex_table_head(caption, 'tbl:gc_mean1', 'l|l|c|c|c|c|c|c|')
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

%% print out a combined info

data_mean = mean(data(:,4:6));

%% save the results in a mat file for later processing

disp('Saving variables...');
save_name = fullfile(log_path,strcat('dfd_gc_results.mat'));

save(save_name, 'data', 'data_name', 'data_sort', 'name_sort', 's_i'); 

%% this section will then start to stitch together the images

generate_images = false;

if(generate_images)

img_space_h = 6;     % number of pixels to put between the images
img_space_v = 6;     % number of pixels to put between the images

data_dir = 'D:\IUPUI\Test_Data\Middlebury_Images_Third';
tb = s_i([1:5,length(data)-4:length(data)]);

% image order: infocus, defocused, groundtruth, result
combined_img = {};

    img_h = 110;
    img_w = 125;

    for idx=1:length(tb)


        %image_files

    %     i_img = imread(fullfile(data_dir, infocus_img_file{tb(idx),1}));
    %     d_img = imread(fullfile(data_dir, defocus_img_file{tb(idx),1}));
    %     gt_img = imread(fullfile(data_dir, gt_img_file{tb(idx),1}));
        i_img = imread(fullfile(data_dir, image_files{tb(idx),1}{1}));
        d_img = imread(fullfile(data_dir, image_files{tb(idx),1}{2}));
        gt_img = imread(fullfile(data_dir, image_files{tb(idx),1}{3}));    
        gt_img = cat(3, gt_img, gt_img, gt_img);

        dm_img = imread(fullfile(image_files{tb(idx),1}{4}));

        % resize the images
        i_img = imresize(i_img,[img_h,img_w]);
        d_img = imresize(d_img,[img_h,img_w]);
        gt_img = imresize(gt_img,[img_h,img_w]);
        dm_img = imresize(dm_img,[img_h,img_w]);

        im_size = size(i_img);

        pad_img = 255*ones([im_size(1) img_space_h im_size(3)]);

        combined_img{idx} = [i_img pad_img d_img pad_img gt_img pad_img dm_img];

    %     combined_img([1:im_size(1), 1:im_size(2), 1:im_size(3)]) = i_img;
    %     combined_img([1:im_size(1), im_size(2)+img_space_h:, 1:im_size(3)]) = i_img;
    %     combined_img([1:im_size(1), 2*1:im_size(2)+img_space_h, 1:im_size(3)]) = i_img;

        %save_dir = 'D:\IUPUI\PhD\Images\dfd_gc';

        save_file = strcat(data_name{tb(idx)}, '_I', num2str(data(tb(idx),1)), '_E', num2str(data(tb(idx),2)-1),'_V',num2str(data(tb(idx),3)),'_combined.png');

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

    save_file = strcat('gc_comb_tb_mb.png');
    imwrite(full_img,fullfile(log_path,save_file));

end



