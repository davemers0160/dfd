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

startpath = 'D:\IUPUI\PhD\Results\k-fold';
[data_file, data_path] = uigetfile(file_filter, 'Select Configuration File', startpath);
if(data_path == 0)
    return;
end

%% Get the data from the file

data_params = parse_input_parameters(fullfile(data_path, data_file));

k_folds = size(data_params,1);

nmae = zeros(k_folds,2);
nrmse = zeros(k_folds,2);
ssim = zeros(k_folds,2);

nmae = zeros(k_folds,1);
nrmse = zeros(k_folds,1);
ssim = zeros(k_folds,1);

for idx=1:k_folds
    nmae(idx,1) = str2double(data_params{idx}{2});
    nrmse(idx,1) = str2double(data_params{idx}{3});
    ssim(idx,1) = str2double(data_params{idx}{4});
    nmae(idx,2) = str2double(data_params{idx}{5});
    nrmse(idx,2) = str2double(data_params{idx}{6});
    ssim(idx,2) = str2double(data_params{idx}{7}); 
end

nmae_stats = zeros(4,2);
nrmse_stats = zeros(4,2);
ssim_stats = zeros(4,2);

nmae_stats(1,1) = min(nmae(:,1));
nmae_stats(2,1) = mean(nmae(:,1));
nmae_stats(3,1) = max(nmae(:,1));
nmae_stats(4,1) = std(nmae(:,1));

nrmse_stats(1,1) = min(nrmse(:,1));
nrmse_stats(2,1) = mean(nrmse(:,1));
nrmse_stats(3,1) = max(nrmse(:,1));
nrmse_stats(4,1) = std(nrmse(:,1));

ssim_stats(1,1) = min(ssim(:,1));
ssim_stats(2,1) = mean(ssim(:,1));
ssim_stats(3,1) = max(ssim(:,1));
ssim_stats(4,1) = std(ssim(:,1));

nmae_stats(1,2) = min(nmae(:,2));
nmae_stats(2,2) = mean(nmae(:,2));
nmae_stats(3,2) = max(nmae(:,2));
nmae_stats(4,2) = std(nmae(:,2));

nrmse_stats(1,2) = min(nrmse(:,2));
nrmse_stats(2,2) = mean(nrmse(:,2));
nrmse_stats(3,2) = max(nrmse(:,2));
nrmse_stats(4,2) = std(nrmse(:,2));

ssim_stats(1,2) = min(ssim(:,2));
ssim_stats(2,2) = mean(ssim(:,2));
ssim_stats(3,2) = max(ssim(:,2));
ssim_stats(4,2) = std(ssim(:,2));

%% plot the training results 
lw = 1.5;
title_str = {'DFD-Net 9-Fold Training Results', 'DFD-Net 9-Fold Test Results'};
plot_save_name = {'_v14a_train', '_v14a_test'};

nrmse_plt_max = ceil(max(nrmse(:))*50)/50;
nmae_plt_max = ceil(max(nmae(:))*50)/50;

    
    fig = figure(plot_num);
    set(gcf,'position',([50,50,1000,700]),'color','w')

    % Plot the NRMSE values
    subplot(3,1,1);
    hold on
    box on
    grid on

    plot([1:k_folds], nrmse(:,1) ,'--b', 'LineWidth', lw);
    plot([1:k_folds], nrmse(:,2) ,'-b', 'LineWidth', lw);
    set(gca, 'fontweight', 'bold', 'FontSize', 13);

    % X-Axis
    xlim([1 k_folds]);
    xticks([1:k_folds]);
    xticklabels([]);

    % Y-Axis
    ylim([0 nrmse_plt_max]);
    yticks([0:0.03:nrmse_plt_max]);
    ytickformat('%1.2f')
    ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

    title(strcat('DFD-Net 9-Fold Preformance Results'),'fontweight','bold','FontSize', 16)

    ax = gca;
    ax.Position = [0.07 0.700 0.915 0.24];

    %--------------------------------------
    % Plot the NMAE values
    subplot(3,1,2);
    hold on
    box on
    grid on

    plot([1:k_folds], nmae(:,1) ,'--b', 'LineWidth', lw);
    plot([1:k_folds], nmae(:,2) ,'-b', 'LineWidth', lw);
    set(gca,'fontweight','bold','FontSize', 13);

    % X-Axis
    xlim([1 k_folds]);
    xticks([1:k_folds]);
    xticklabels([]);

    % Y-Axis
    ylim([0 nmae_plt_max]);
    yticks([0:0.02:nmae_plt_max]);
    ytickformat('%1.2f');
    ylabel('NMAE','fontweight','bold','FontSize', 13);

    %title(strcat('DFD-Net Average Training Performance Results'),'fontweight','bold','FontSize', 16)

    ax = gca;
    ax.Position = [0.07 0.425 0.915 0.24];

    %--------------------------------------
    % Plot the SSIM values
    subplot(3,1,3);
    hold on
    box on
    grid on

    plot([1:k_folds], ssim(:,1) ,'--b', 'LineWidth', lw);
    plot([1:k_folds], ssim(:,2) ,'-b', 'LineWidth', lw);
    set(gca,'fontweight','bold','FontSize', 13);

    % X-Axis
    xlim([1 k_folds]);
    xticks([1:k_folds]);
    %xtickangle(90);
    %xtickformat();
    xlabel(strcat('9-Fold Test Number'), 'fontweight', 'bold', 'FontSize', 13);

    % Y-Axis
    ylim([0.7 1]);
    yticks([0.7:0.1:1]);
    ytickformat('%1.2f');
    ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);
    legend('Training Results', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')

    ax = gca;
    ax.Position = [0.07 0.145 0.915 0.24];

    print(plot_num, '-dpng', fullfile(data_path,strcat('dfd_dnn_k-fold_v14a_combined','.png')));

    plot_num = plot_num + 1;  

%% create a table for latex
fprintf('\n\n');

caption = {'DFD-Net Average Training Performance Results', 'DFD-Net Average Test Performance Results'};
row_name = {'Minimum', 'Mean', 'Maximum', 'Std Deviation'};

write_latex_table_head('DFD-Net Average 9-Fold Cross Validation Performance Results', 'tbl:dnn_k_fold1', 'l|c|c|c|c|c|c|')
fprintf('\\cline{2-7} \n');
fprintf('\\multicolumn{1}{c|}{} & \\multicolumn{3}{c|}{\\textbf{Training}} & \\multicolumn{3}{c|}{\\textbf{Testing}} \\\\ \\cline{2-7} \n');
fprintf('\\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{\\textbf{NRMSE}} & \\multicolumn{1}{c|}{\\textbf{NMAE}} & \\multicolumn{1}{c|}{\\textbf{SSIM}} & \\multicolumn{1}{c|}{\\textbf{NRMSE}} & \\multicolumn{1}{c|}{\\textbf{NMAE}} & \\multicolumn{1}{c|}{\\textbf{SSIM}} \\\\ \\hline \n');

% print out the actual values
for idx=1:k_folds
    fprintf('\\multicolumn{1}{|l|}{\\textbf{Test %d}} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', idx, nrmse(idx,1), nmae(idx,1), ssim(idx,1), nrmse(idx,2), nmae(idx,2), ssim(idx,2));   
end

% print out the stat values
for idx=1:4
    fprintf('\\multicolumn{1}{|l|}{\\textbf{%s}} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', row_name{idx}, nrmse_stats(idx,1), nmae_stats(idx,1), ssim_stats(idx,1), nrmse_stats(idx,2), nmae_stats(idx,2), ssim_stats(idx,2));   
end
write_latex_table_end;

fprintf('\n\n');


