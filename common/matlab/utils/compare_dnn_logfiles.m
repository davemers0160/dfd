format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

%% setup some variables
plot_num = 1;
sep = '------------------------------------------------------------------';
file_filter = {'*.txt','Text Files';'*.*','All Files' };
lw = 1.00;
lw_stem = 1.0;

%% get the log file to parse
if(ispc)
    startpath = 'D:\IUPUI\PhD\Results';
else
    startpath = '/home/owner/DfD';
end

[log_file1, log_path1] = uigetfile(file_filter, 'Select First DFD-Net Log File', startpath);
if(log_path1 == 0)
    return;
end

[log_file2, log_path2] = uigetfile(file_filter, 'Select Second DFD-Net Log File', log_path1);
if(log_path2 == 0)
    return;
end

save_path = 'D:\IUPUI\PhD\Images\dfd_dnn';

commandwindow;

%% Start reading the log file
[log_version{1}, data{1}, final_traning{1}, final_test{1}] = read_dnn_logfile(fullfile(log_path1,log_file1));

[log_version{2}, data{2}, final_traning{2}, final_test{2}] = read_dnn_logfile(fullfile(log_path2,log_file2));

[lr{1}, lr_idx{1}] = unique([data{1}.lr]','stable');
[lr{2}, lr_idx{2}] = unique([data{2}.lr]','stable');

steps{1} = [data{1}.step]';
steps{2} = [data{2}.step]';

min_steps = min([min(steps{1}), min(steps{2})]);
max_steps = max([max(steps{1}), max(steps{2})]);

legend_names = {'Basic Method Training', 'Basic Method Test', 'Basic Method Convergence Step', 'Enhanced Method Training', 'Enhanced Method Test', 'Enhanced Method Convergence Step'};

%% plot the NMAE Comparisons

max_nmae = max([max([data{1}.tr_nmae]),max([data{1}.te_nmae]),max([data{2}.tr_nmae]),max([data{2}.te_nmae])]);
max_nmae = (ceil(max_nmae*20)/20);

figure(plot_num)
set(gcf,'position',([100,100,1300,600]),'color','w')
hold on
box on
grid on
% plot the first set
p1 = plot(steps{1}, [data{1}.tr_nmae]','-b');
p2 = plot(steps{1}, [data{1}.te_nmae]','-g');
p3 = stem(steps{1}(lr_idx{1}(2:end)), max_nmae*ones(1,numel(lr_idx{1}(2:end))),'LineStyle','-','Color','r','Marker','none');
%plot the second set
p4 = plot(steps{2}, [data{2}.tr_nmae]','--b');
p5 = plot(steps{2}, [data{2}.te_nmae]','--g');
p6 = stem(steps{2}(lr_idx{2}(2:end)), max_nmae*ones(1,numel(lr_idx{2}(2:end))),'LineStyle','--','Color','r','Marker','none');

% X-Axis
set(gca, 'fontweight','bold')
xlim([min_steps max_steps]);
%xticks([1:length(E_name{1})]);
%xticklabels(E_name{1});
xlabel('Training Steps', 'fontweight','bold');

% Y-Axis
ylim([0 max_nmae]);
ylabel('NMAE', 'fontweight','bold');
ytickformat('%1.2f')

title('Training and Testing NMAE Results', 'fontweight','bold');
legend([p1 p4 p2 p5 p3 p6], legend_names, 'location', 'southoutside', 'orientation','horizontal');

ax = gca;
ax.Position = [0.05 0.15 0.93 0.80];


plot_num = plot_num + 1;

%% plot the NRMSE comparisons

max_nrmse = max([max([data{1}.tr_nrmse]),max([data{1}.te_nrmse]),max([data{2}.tr_nrmse]),max([data{2}.te_nrmse])]);
max_nrmse = (ceil(max_nrmse*20)/20);

figure(plot_num)
set(gcf,'position',([100,100,1300,600]),'color','w')
hold on
box on
grid on
% plot the first set
p1 = plot(steps{1}, [data{1}.tr_nrmse]','-b');
p2 = plot(steps{1}, [data{1}.te_nrmse]','-g');
p3 = stem(steps{1}(lr_idx{1}(2:end)), max_nrmse*ones(1,numel(lr_idx{1}(2:end))),'LineStyle','-','Color','r','Marker','none');
%plot the second set
p4 = plot(steps{2}, [data{2}.tr_nrmse]','--b');
p5 = plot(steps{2}, [data{2}.te_nrmse]','--g');
p6 = stem(steps{2}(lr_idx{2}(2:end)), max_nrmse*ones(1,numel(lr_idx{2}(2:end))),'LineStyle','--','Color','r','Marker','none') ;

% X-Axis
set(gca, 'fontweight','bold')
xlim([min_steps max_steps]);
%xticks([1:length(E_name{1})]);
%xticklabels(E_name{1});
xlabel('Training Steps', 'fontweight','bold');

% Y-Axis
ylim([0 max_nrmse]);
ylabel('NRMSE', 'fontweight','bold');
ytickformat('%1.2f')

title('Training and Testing NRMSE Results', 'fontweight','bold');
legend([p1 p4 p2 p5 p3 p6], legend_names, 'location', 'southoutside', 'orientation','horizontal');

ax = gca;
ax.Position = [0.05 0.15 0.93 0.80];

plot_num = plot_num + 1;


%% plot the SSIM comparisons

min_ssim = 0;
max_ssim = 1;


figure(plot_num)
set(gcf,'position',([100,100,1300,600]),'color','w')

hold on
box on
grid on
% plot the first set
p1 = plot(steps{1}, [data{1}.tr_ssim]','-b');
p2 = plot(steps{1}, [data{1}.te_ssim]','-g');
p3 = stem(steps{1}(lr_idx{1}(2:end)), max_ssim*ones(1,numel(lr_idx{1}(2:end))),'LineStyle','-','Color','r','Marker','none');
%plot the second set
p4 = plot(steps{2}, [data{2}.tr_ssim]','--b');
p5 = plot(steps{2}, [data{2}.te_ssim]','--g');
p6 = stem(steps{2}(lr_idx{2}(2:end)), max_ssim*ones(1,numel(lr_idx{2}(2:end))),'LineStyle','--','Color','r','Marker','none') ;

% X-Axis
set(gca, 'fontweight','bold')
xlim([min_steps max_steps]);
%xticks([1:length(E_name{1})]);
%xticklabels(E_name{1});
xlabel('Training Steps', 'fontweight','bold');

% Y-Axis
ylim([min_ssim max_ssim]);
ylabel('SSIM', 'fontweight','bold');
ytickformat('%1.2f')

title('Training and Testing SSIM Results', 'fontweight','bold');

legend([p1 p4 p2 p5 p3 p6], legend_names, 'location', 'southoutside', 'orientation','horizontal');

ax = gca;
ax.Position = [0.05 0.15 0.93 0.80];

plot_num = plot_num + 1;

%% Put all of the graphs on one plot

% NRMSE
figure(plot_num)
set(gcf,'position',([50,50,1300,700]),'color','w')
subplot(3,1,1);
hold on
box on
grid on
% plot the first set
p1 = plot(steps{1}, [data{1}.tr_nrmse]','-b', 'LineWidth', lw);
p2 = plot(steps{1}, [data{1}.te_nrmse]','-g', 'LineWidth', lw);
%p3 = stem(steps{1}(lr_idx{1}(2:end)), max_nrmse*ones(1,numel(lr_idx{1}(2:end))),'LineStyle','-','Color','r','Marker','none', 'LineWidth', lw_stem);
p3 = stem(steps{1}(end), max_nrmse, 'LineStyle','-', 'Color','r', 'Marker','none', 'LineWidth', lw_stem);
%plot the second set
p4 = plot(steps{2}, [data{2}.tr_nrmse]','--b', 'LineWidth', lw);
p5 = plot(steps{2}, [data{2}.te_nrmse]','--g', 'LineWidth', lw);
%p6 = stem(steps{2}(lr_idx{2}(2:end)), max_nrmse*ones(1,numel(lr_idx{2}(2:end))),'LineStyle','--','Color','r','Marker','none', 'LineWidth', lw_stem) ;
p6 = stem(steps{2}(end), max_nrmse, 'LineStyle','--', 'Color','r', 'Marker','none', 'LineWidth', lw_stem) ;

% X-Axis
set(gca, 'fontweight','bold', 'fontsize',13);
xlim([min_steps max_steps]);
xticklabels([]);
%xticks([1:length(E_name{1})]);
%xlabel('Training Steps', 'fontweight','bold');

% Y-Axis
ylim([0 max_nrmse]);
yticks([0:0.1:0.5]);
ylabel('NRMSE', 'fontweight','bold');
ytickformat('%1.2f')

title('Training Method Comparison', 'fontweight','bold', 'fontsize',16);

ax = gca;
ax.Position = [0.06 0.71 0.93 0.24];

% NMAE
subplot(3,1,2)
hold on
box on
grid on
% plot the first set
p1 = plot(steps{1}, [data{1}.tr_nmae]','-b', 'LineWidth', lw);
p2 = plot(steps{1}, [data{1}.te_nmae]','-g', 'LineWidth', lw);
%p3 = stem(steps{1}(lr_idx{1}(2:end)), max_nmae*ones(1,numel(lr_idx{1}(2:end))),'LineStyle','-','Color','r','Marker','none', 'LineWidth', lw_stem);
p3 = stem(steps{1}(end), max_nmae, 'LineStyle','-', 'Color','r', 'Marker','none', 'LineWidth', lw_stem);
%plot the second set
p4 = plot(steps{2}, [data{2}.tr_nmae]','--b', 'LineWidth', lw);
p5 = plot(steps{2}, [data{2}.te_nmae]','--g', 'LineWidth', lw);
%p6 = stem(steps{2}(lr_idx{2}(2:end)), max_nmae*ones(1,numel(lr_idx{2}(2:end))),'LineStyle','--','Color','r','Marker','none', 'LineWidth', lw_stem);
p6 = stem(steps{2}(end), max_nmae, 'LineStyle','--', 'Color','r', 'Marker','none', 'LineWidth', lw_stem);

% X-Axis
set(gca, 'fontweight','bold', 'fontsize',13);
xlim([min_steps max_steps]);
xticklabels([]);
%xticks([1:length(E_name{1})]);
%xlabel('Training Steps', 'fontweight','bold');

% Y-Axis
ylim([0 max_nmae]);
yticks([0:0.1:0.5]);
ylabel('NMAE', 'fontweight','bold');
ytickformat('%1.2f')

%title('Training and Testing NMAE Results', 'fontweight','bold');

ax = gca;
ax.Position = [0.06 0.435 0.93 0.24];

% plot the SSIM comparisons
min_ssim = 0;
max_ssim = 1;

% figure(plot_num)
% set(gcf,'position',([100,100,1400,600]),'color','w')
subplot(3,1,3);
hold on
box on
grid on
% plot the first set
p1 = plot(steps{1}, [data{1}.tr_ssim]','-b', 'LineWidth', lw);
p2 = plot(steps{1}, [data{1}.te_ssim]','-g', 'LineWidth', lw);
%p3 = stem(steps{1}(lr_idx{1}(2:end)), max_ssim*ones(1,numel(lr_idx{1}(2:end))),'LineStyle','-','Color','r','Marker','none', 'LineWidth', lw_stem);
p3 = stem(steps{1}(end), max_ssim, 'LineStyle','-', 'Color','r', 'Marker','none', 'LineWidth', lw_stem);
%plot the second set
p4 = plot(steps{2}, [data{2}.tr_ssim]','--b', 'LineWidth', lw);
p5 = plot(steps{2}, [data{2}.te_ssim]','--g', 'LineWidth', lw);
%p6 = stem(steps{2}(lr_idx{2}(2:end)), max_ssim*ones(1,numel(lr_idx{2}(2:end))),'LineStyle','--','Color','r','Marker','none', 'LineWidth', lw_stem) ;
p6 = stem(steps{2}(end), max_ssim, 'LineStyle','--', 'Color','r', 'Marker','none', 'LineWidth', lw_stem) ;

% X-Axis
set(gca, 'fontweight','bold','fontsize',13)
xlim([min_steps max_steps]);
%xticks([1:length(E_name{1})]);
%xticklabels(E_name{1});
xlabel('Training Steps', 'fontweight','bold');

% Y-Axis
ylim([min_ssim max_ssim]);
yticks([0:0.25:1.0]);
ylabel('SSIM', 'fontweight','bold');
ytickformat('%1.2f')

%title('Training and Testing SSIM Results', 'fontweight','bold');
lgd = legend([p1 p2 p3 p4 p5 p6], legend_names, 'location', 'southoutside', 'orientation','horizontal');
lgd.NumColumns = 3;


ax = gca;
ax.Position = [0.06 0.16 0.93 0.24];

print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_net_v14a_tng_patch_comp.png')));

plot_num = plot_num + 1;


