format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

plot_num = 1;

%% get the two mat files

% expect that the histograms are in the raw form for each mat file

start_path = 'D:\IUPUI\PhD\';
%save_path = 'D:\IUPUI\PhD\Images\datasets';

file_filter = {'*.mat','MAT-files (*.mat)'; '*.*','All Files' };

[d1_file, d1_path] = uigetfile(file_filter,'Select First Distribution Mat File', start_path);
if(d1_path == 0)
    return;
end

[d2_file, d2_path] = uigetfile(file_filter,'Select Second Distribution Mat File', start_path);
if(d2_path == 0)
    return;
end

% get the save file name
file_filter = {'*.png','Image Files';'*.*','All Files' };
[save_name, save_path] = uiputfile(file_filter, 'Enter save file', start_path);

if(save_path == 0)
    save_image = false;
else
    save_image = true;
end


commandwindow;

%% load in the data


d1 = load(fullfile(d1_path, d1_file));
d2 = load(fullfile(d2_path, d2_file));

% Only use the training hist
sum_d1_hist = sum(d1.train_hist,1);
sum_d2_hist = sum(d2.train_hist,1);


% calculate the ratio between the two
d1_d2_ratio = sum_d2_hist./sum_d1_hist;
nn = ~isnan(d1_d2_ratio);

fprintf('Min Ratio:     %2.2f\n', max(d1_d2_ratio(nn)));
fprintf('Average Ratio: %2.2f\n', mean(d1_d2_ratio(nn)));
fprintf('Max Ratio:     %2.2f\n', min(d1_d2_ratio(nn)));

%% Plot the overlay

sum_d1_train_hist = sum(d1.train_hist_sum(:));
sum_d2_train_hist = sum(d2.train_hist_sum(:));

min_bin = 0;
max_bin = numel(sum_d1_hist);

x_lim = [0,255];

hist_bin_step = 1;
hist_bins = min_bin:hist_bin_step:max_bin;

plot_step = 5;
%combined_hist = [sum_d1_hist/sum(sum_d1_hist); sum_d2_hist/sum(sum_d2_hist)]';
%combined_hist = [sum_d1_hist/sum_d1_train_hist; sum_d2_hist/sum_d2_train_hist]';
combined_hist = [sum_d1_hist; sum_d2_hist]';

x = [min_bin:hist_bin_step:(max_bin-1)];
y_m = ceil(log10(max(combined_hist(:))));
y_max = 10^y_m;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
b = bar(x, combined_hist);
set(gca,'fontweight','bold','FontSize',13, 'yscale', 'log');

% X-Axis
xlim(x_lim);
%xticks([(x_lim(1)+2):plot_step:(x_lim(2)-2)]);
xticks([0:plot_step:max_bin]);
xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ylabel('Depth Map Value Count','fontweight','bold')
ylim([1 y_max]);
%ytickformat('%1.2f');

b(1).FaceColor = 'b';
b(2).FaceColor = 'r';

title('Depth Map Distribution Comparison', 'fontweight','bold','FontSize',16);
lgd = legend(strcat(d1.name,32,'Distribution'),strcat(d2.name,32,'Distribution'), 'location','southoutside', 'orientation', 'horizontal');
ax = gca;
ax.Position = [0.07 0.17 0.91 0.77];


%save_name = strcat(d1.name,'_',d2.name,'_hist_comp.png');
if(save_image)
    print(plot_num, '-dpng', fullfile(save_path, save_name));
end
plot_num = plot_num + 1;




