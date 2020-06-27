format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% select the data file

file_filter = {'*.txt','Text Files';'*.*','All Files' };

startpath = 'D:\IUPUI\PhD\Results\color_v_patch';
[data_file, data_path] = uigetfile(file_filter, 'Select the Training Patch Data File', startpath);
if(data_path == 0)
    return;
end

startpath = 'D:\IUPUI\PhD\Results\dfd_dnn_trials';
[multi_results_file, results_path] = uigetfile(file_filter, 'Select Multiple Trials Results File', startpath);
if(results_path == 0)
    return;
end

commandwindow;

save_path = 'D:\IUPUI\PhD\Images\dfd_dnn';

%% Get the data from the training patch file

data_params = parse_input_parameters(fullfile(data_path, data_file));

num_trials = size(data_params,1);

nmae = zeros(num_trials,2);
nrmse = zeros(num_trials,2);
ssim = zeros(num_trials,2);
patch_label = cell(num_trials,1);

for idx=1:num_trials
    patch_label{idx,1} = data_params{idx}{1};
    nmae(idx,1) = str2double(data_params{idx}{2});
    nrmse(idx,1) = str2double(data_params{idx}{3});
    ssim(idx,1) = str2double(data_params{idx}{4});
    nmae(idx,2) = str2double(data_params{idx}{5});
    nrmse(idx,2) = str2double(data_params{idx}{6});
    ssim(idx,2) = str2double(data_params{idx}{7});
end


%% get the data from the multiple training events file

data_params = parse_input_parameters(fullfile(results_path, multi_results_file));

num_trials2 = size(data_params,1);

mt_nmae = zeros(num_trials2,2);
mt_nrmse = zeros(num_trials2,2);
mt_ssim = zeros(num_trials2,2);

mt_nmae = zeros(num_trials2,1);
mt_nrmse = zeros(num_trials2,1);
mt_ssim = zeros(num_trials2,1);

for idx=1:num_trials2
    mt_nmae(idx,1) = str2double(data_params{idx}{1});
    mt_nrmse(idx,1) = str2double(data_params{idx}{2});
    mt_ssim(idx,1) = str2double(data_params{idx}{3});
    mt_nmae(idx,2) = str2double(data_params{idx}{4});
    mt_nrmse(idx,2) = str2double(data_params{idx}{5});
    mt_ssim(idx,2) = str2double(data_params{idx}{6}); 
end

nmae_stats = zeros(4,2);
nrmse_stats = zeros(4,2);
ssim_stats = zeros(4,2);

nmae_stats(1,1) = min(mt_nmae(:,1));
nmae_stats(2,1) = mean(mt_nmae(:,1));
nmae_stats(3,1) = max(mt_nmae(:,1));
nmae_stats(4,1) = std(mt_nmae(:,1));

nrmse_stats(1,1) = min(mt_nrmse(:,1));
nrmse_stats(2,1) = mean(mt_nrmse(:,1));
nrmse_stats(3,1) = max(mt_nrmse(:,1));
nrmse_stats(4,1) = std(mt_nrmse(:,1));

ssim_stats(1,1) = min(mt_ssim(:,1));
ssim_stats(2,1) = mean(mt_ssim(:,1));
ssim_stats(3,1) = max(mt_ssim(:,1));
ssim_stats(4,1) = std(mt_ssim(:,1));

nmae_stats(1,2) = min(mt_nmae(:,2));
nmae_stats(2,2) = mean(mt_nmae(:,2));
nmae_stats(3,2) = max(mt_nmae(:,2));
nmae_stats(4,2) = std(mt_nmae(:,2));

nrmse_stats(1,2) = min(mt_nrmse(:,2));
nrmse_stats(2,2) = mean(mt_nrmse(:,2));
nrmse_stats(3,2) = max(mt_nrmse(:,2));
nrmse_stats(4,2) = std(mt_nrmse(:,2));

ssim_stats(1,2) = min(mt_ssim(:,2));
ssim_stats(2,2) = mean(mt_ssim(:,2));
ssim_stats(3,2) = max(mt_ssim(:,2));
ssim_stats(4,2) = std(mt_ssim(:,2));

nrmse_plt_max = ceil(max(mt_nrmse(:))*25)/25;
nmae_plt_max = ceil(max(mt_nmae(:))*34)/34;

lw = 1.5;
sig_m = 2;

%% combined trainging and testing plot

if(false)
    
    fig = figure(plot_num);
    set(gcf,'position',([50,50,1000,700]),'color','w')

    % Plot the NRMSE values
    subplot(3,1,1);
    hold on
    box on
    grid on

    plot([1:num_trials], nrmse(:,1) ,'--b', 'LineWidth', lw);
    plot([1:num_trials], nrmse(:,2) ,'-b', 'LineWidth', lw);
    set(gca, 'fontweight', 'bold', 'FontSize', 13);

    % X-Axis
    xlim([1 num_trials]);
    xticks([1:num_trials]);
    xticklabels([]);

    % Y-Axis
    ylim([0 nrmse_plt_max]);
    yticks([0:0.1:nrmse_plt_max]);
    ytickformat('%1.2f')
    ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

    title('Training Patch Size Performance Comparison','fontweight','bold','FontSize', 16)

    ax = gca;
    ax.Position = [0.07 0.710 0.92 0.23];

    %--------------------------------------
    % Plot the NMAE values
    subplot(3,1,2);
    hold on
    box on
    grid on

    plot([1:num_trials], nmae(:,1) ,'--b', 'LineWidth', lw);
    plot([1:num_trials], nmae(:,2) ,'-b', 'LineWidth', lw);
    set(gca,'fontweight','bold','FontSize', 13);

    % X-Axis
    xlim([1 num_trials]);
    xticks([1:num_trials]);
    xticklabels([]);

    % Y-Axis
    ylim([0 nmae_plt_max]);
    yticks([0:0.1:nmae_plt_max]);
    ytickformat('%1.2f');
    ylabel('NMAE','fontweight','bold','FontSize', 13);

    %title(strcat('DFD-Net Average Training Performance Results'),'fontweight','bold','FontSize', 16)

    ax = gca;
    ax.Position = [0.07 0.445 0.92 0.23];

    %--------------------------------------
    % Plot the SSIM values
    subplot(3,1,3);
    hold on
    box on
    grid on

    plot([1:num_trials], ssim(:,1) ,'--b', 'LineWidth', lw);
    plot([1:num_trials], ssim(:,2) ,'-b', 'LineWidth', lw);
    set(gca,'fontweight','bold','FontSize', 13);

    % X-Axis
    xlim([1 num_trials]);
    xticks([1:num_trials]);
    xtickangle(90);
    %xtickformat();
    xticklabels(patch_size);
    xlabel(strcat('Training Patch Size'), 'fontweight', 'bold', 'FontSize', 13);

    % Y-Axis
    ylim([0.2 1]);
    yticks([0:0.2:1]);
    ytickformat('%1.2f');
    ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);
    legend('Training Results', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')

    ax = gca;
    ax.Position = [0.07 0.175 0.92 0.23];

    print(plot_num, '-dpng', fullfile(data_path,strcat('dfd_dnn_combine_patch_size.png')));

    plot_num = plot_num + 1;  
    
end

%% Plot the path size data


% get the lower limits
limits_test(1,1) = nmae_stats(2,2) - sig_m*nmae_stats(4,2);
limits_test(2,1) = nrmse_stats(2,2) - sig_m*nrmse_stats(4,2);
limits_test(3,1) = ssim_stats(2,2) - sig_m*ssim_stats(4,2);

% get the upper limits
limits_test(1,2) = nmae_stats(2,2) + sig_m*nmae_stats(4,2);
limits_test(2,2) = nrmse_stats(2,2) + sig_m*nrmse_stats(4,2);
limits_test(3,2) = ssim_stats(2,2) + sig_m*ssim_stats(4,2);


figure(plot_num)
set(gcf,'position',([50,50,1200,700]),'color','w')
% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
area([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], 'BaseValue', limits_test(2,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth', 1);

%plot the mean line - maybe
plot([-1,num_trials+2],[nrmse_stats(2,2), nrmse_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], nrmse(:,2), 18,'filled','k');
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
xticklabels([]);

% Y-Axis
%nrmse_plt_max = ceil(limits_test(2,2)*80)/80;
%nrmse_plt_min = floor(limits_test(2,1)*80)/80;
nrmse_plt_max = 0.09;
nrmse_plt_min = 0.05;
ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:0.01:nrmse_plt_max]);
ytickformat('%1.3f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

title('Training Patch Size Performance Comparison','fontweight','bold','FontSize', 16)

ax = gca;
ax.Position = [0.070 0.710 0.92 0.23];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
area([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], 'BaseValue', limits_test(1,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1);

%plot the mean line - maybe
plot([-1,num_trials+2],[nmae_stats(2,2), nmae_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], nmae(:,2), 18,'filled','k');
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
xticklabels([]);

% Y-Axis
%nmae_plt_max = ceil(limits_test(1,2)*200)/200;
%nmae_plt_min = floor(limits_test(1,1)*200)/200;
nmae_plt_max = 0.038;
nmae_plt_min = 0.014;
ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:0.004:nmae_plt_max]);
ytickformat('%1.3f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
ax.Position = [0.070 0.445 0.92 0.23];
%--------------------------------------
subplot(3,1,3);
hold on
box on
grid on

%plot the area
p1 = area([-1,num_trials+2],[limits_test(3,2), limits_test(3,2)], 'BaseValue', limits_test(3,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1);

%plot the mean line - maybe
p2 = plot([-1,num_trials+2],[ssim_stats(2,2), ssim_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(3,2), limits_test(3,2)], '-b', 'LineWidth', 1);

% plot the test points
p3 = scatter([1:num_trials], ssim(:,2), 18,'filled','k');
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
% xlabel_name = {''};
% for idx=1:num_trials
%     xlabel_name{idx+1,1} = num2str(idx);
% end
xlim([0 num_trials+1]);
xticks([1:num_trials]);
xticklabels(patch_label);
xtickangle(90);

xlabel(strcat('Patch Size'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
%ssim_plt_max = ceil(limits_test(3,2)*80)/80;
%ssim_plt_min = floor(limits_test(3,1)*80)/80;
ssim_plt_max = 0.925;
ssim_plt_min = 0.80;
ylim([ssim_plt_min ssim_plt_max]);
yticks([ssim_plt_min:0.025:ssim_plt_max]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
ytickformat('%1.3f')
ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

legend([p1, p2, p3], strcat(32,'{\pm}',num2str(sig_m),'{\sigma} Metric Bounds'), 'Metric Mean', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.070 0.175 0.92 0.23];


print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_dnn_patch_size_results.png')));

plot_num = plot_num + 1;  


