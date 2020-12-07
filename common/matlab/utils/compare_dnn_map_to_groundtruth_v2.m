format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% select the data file with the images

file_filter = {'*.txt','Text Files';'*.*','All Files' };
startpath = 'D:\Projects\';

[gt_input_file, gt_input_path] = uigetfile(file_filter, 'Select Train/Test Input File', startpath);
if(gt_input_path == 0)
    return;
end

[dm_input_file, dm_input_path] = uigetfile(file_filter, 'Select Depthmap Results File', startpath);
if(dm_input_path == 0)
    return;
end

commandwindow;

%% check to see if the inputs are normal or grouped

gt_params = parse_input_parameters(fullfile(gt_input_path, gt_input_file));

% get the directory for the data
gt_data_directory = gt_params{1}{1};
gt_params(1) = [];


dm_params = parse_input_parameters(fullfile(dm_input_path, dm_input_file));

if(size(gt_params,1) ~= size(dm_params,1))
    num = min(size(gt_params,1), size(dm_params,1));
else
    num = size(gt_params,1);
end

%% get the maximum depthmap value
prompt = {'Enter maximum depthmap value:'};
dlgtitle = '';
definput = {'20'};
answer = inputdlg(prompt,dlgtitle,[1 40],definput);

max_depthmap_value = str2double(answer);

%% run through each depth map and get the results

combined_hist = zeros(max_depthmap_value+1, 2);
cm = zeros(max_depthmap_value+1, max_depthmap_value+1);

for idx=1:num

    gt_img = imread(fullfile(gt_data_directory, gt_params{idx}{3}));
    
    dnn_img = imread(dm_params{idx}{1});
    dnn_img_size = size(dnn_img);
    
    [gt_error_value, tmp_combined_hist, tmp_cm] = get_depthmap_errors(gt_img, dnn_img, max_depthmap_value);

    combined_hist = combined_hist + tmp_combined_hist;
    cm = cm  + tmp_cm;

end


sum_cm = sum(cm,2);
cm_diag = diag(cm);

cm_correct = cm_diag./sum_cm;


%% try some plotting

figure(plot_num)
set(gcf,'position',([100,100,1600,800]),'color','w')
cm_plot = confusionchart(cm, [0:1:max_depthmap_value], 'RowSummary','row-normalized','ColumnSummary','column-normalized');
%set(gca,'fontweight','bold','FontSize',13);

cm_plot.FontSize = 11;

xlabel('Predicted Depthmap Values');
ylabel('Actual Depthmap Values');

ax = gca;
ax.Position = [0.04 0.07 0.94 0.9];

print(plot_num, '-dpng', fullfile(dm_input_path, 'depth_map_results_cm.png'));
savefig(gcf, fullfile(dm_input_path, 'depth_map_results_cm.fig'), 'compact');

plot_num = plot_num + 1;

%% plot the percent correct classification

x_lim = [-1,max_depthmap_value+1];
x = 0:1:max_depthmap_value;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
bar(x, (1-cm_correct)*100,'b')
set(gca,'fontweight','bold','FontSize',13);
grid on
box on

% X-Axis
xlim(x_lim);
xticks([0:1:max_depthmap_value]);
xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ylim([0 100]);
yticks([0:10:100]);
ylabel('Depth Map Error Percentage','fontweight','bold')

ax = gca;
ax.Position = [0.065 0.1 0.91 0.86];

print(plot_num, '-dpng', fullfile(dm_input_path, 'depth_map_results_error.png'));
savefig(gcf, fullfile(dm_input_path, 'depth_map_results_error.fig'), 'compact');

plot_num = plot_num + 1;

%% histogram plot

min_bin = 0;
max_bin = size(combined_hist,1);

plot_step = 1;
hist_bin_step = 1;
hist_bins = min_bin:hist_bin_step:max_bin;

x_lim = [min_bin-1,max_depthmap_value+1];
x = [min_bin:hist_bin_step:(max_bin-1)];

y_m = ceil(log10(max(combined_hist(:))));
%y_max = 10^y_m;
y_max = 1.02*max(combined_hist(:))+10;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
b = bar(x, combined_hist);
%set(gca,'fontweight','bold','FontSize',13, 'yscale', 'log');
set(gca,'fontweight','bold','FontSize',13);

% X-Axis
xlim(x_lim);
%xticks([(x_lim(1)+2):plot_step:(x_lim(2)-2)]);
xticks([0:plot_step:max_bin-1]);
xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ylabel('Depth Map Value Count','fontweight','bold')
ylim([1 y_max]);
%ytickformat('%1.2f');

b(1).FaceColor = 'b';
b(2).FaceColor = 'r';

%title('Depth Map Distribution Comparison', 'fontweight','bold','FontSize',16);
lgd = legend('Ground Truth', 'DNN Depth Map', 'location','southoutside', 'orientation', 'horizontal');
ax = gca;
ax.Position = [0.06 0.17 0.92 0.78];

plot_num = plot_num + 1;





