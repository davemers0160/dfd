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

startpath = 'D:\IUPUI\PhD\Results\dfd_dnn_trials';
[data_file, data_path] = uigetfile(file_filter, 'Select Configuration File', startpath);
if(data_path == 0)
    return;
end

run_test_data = true;
startpath = 'D:\IUPUI\PhD\Results\dfd_dnn\dnn_reduction';
test_data_path = uigetdir(startpath, 'Select Multi Trial Results Folder');
if(test_data_path == 0)
    run_test_data = false;
end


save_results = false;

startpath = 'D:\IUPUI\PhD\';
save_path = uigetdir(startpath,'Select Save Folder');

if(save_path ~= 0)
    save_results = true;
end

commandwindow;

%% Get the data from the file

dfd_dnn_index = [1,2,3,4,5,6];
[data_name, num_trials, nmae, nrmse, ssim, nmae_stats, nrmse_stats, ssim_stats] = get_dfd_results(fullfile(data_path, data_file), dfd_dnn_index);

lw = 1.5;
sig_m = 2;
test_sig_m = 2;

fprintf('\n\nName:, %s\n', data_name);
fprintf(' , NMAE, NRMSE, SSIM\n');
fprintf('Min:, %2.4f, %2.4f, %2.4f\n', nmae_stats(1,2), nrmse_stats(1,2), ssim_stats(1,2));
fprintf('Mean:, %2.4f, %2.4f, %2.4f\n', nmae_stats(2,2), nrmse_stats(2,2), ssim_stats(2,2));
fprintf('Max:, %2.4f, %2.4f, %2.4f\n', nmae_stats(3,2), nrmse_stats(3,2), ssim_stats(3,2));
fprintf('Std:, %2.4f, %2.4f, %2.4f\n', nmae_stats(4,2), nrmse_stats(4,2), ssim_stats(4,2));

combined_test = cat(2, nmae(:,2), nrmse(:,2), ssim(:,2));
[comb_test_sort, comb_index] = sortrows(combined_test, [2,1], 'ascend');


%% Another way of showing the same info
% Plot the 2 sigma area along with the results from each of the
% training/test results

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
%set(gcf,'position',([50,50,800,500]),'color','w')

% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(2,2); num_trials+2,limits_test(2,2); num_trials+2,limits_test(2,1); -1,limits_test(2,1)];
patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);

%plot the mean line - maybe
plot([-1,num_trials+2],[nrmse_stats(2,2), nrmse_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], nrmse(:,2),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xticklabels([]);

% Y-Axis
p = 100;
nrmse_plt_max = max(ceil(limits_test(2,2)*p)/p, ceil(nrmse_stats(3,2)*p)/p);
nrmse_plt_min = min(floor(limits_test(2,1)*p)/p, floor(nrmse_stats(1,2)*p)/p);
%dfd
% nrmse_plt_max = 0.08;
% nrmse_plt_min = 0.05;
%dfd rw4
% nrmse_plt_max = 0.075;
% nrmse_plt_min = 0.060;
ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:.01:nrmse_plt_max]);
ytickformat('%1.3f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

%title('DfD-Net Real World - Test Results Distribution','fontweight','bold','FontSize', 16)
title(data_name,'fontweight','bold','FontSize', 16)

ax = gca;
%ax.XAxis.MinorTickValues = 1:1:num_trials;
ax.Position = [0.078 0.700 0.905 0.24];
%ax.Position = [0.10 0.710 0.88 0.23];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(1,2); num_trials+2,limits_test(1,2); num_trials+2,limits_test(1,1); -1,limits_test(1,1)];
patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);

%plot the mean line - maybe
plot([-1,num_trials+2],[nmae_stats(2,2), nmae_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], nmae(:,2),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xticklabels([]);

% Y-Axis
p = 500;
nmae_plt_max = max(ceil(limits_test(1,2)*p)/p, ceil(nmae_stats(3,2)*p)/p);
nmae_plt_min = min(floor(limits_test(1,1)*p)/p, floor(nmae_stats(1,2)*p)/p);
% nmae_plt_max = 0.021;
% nmae_plt_min = 0.015;
ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:0.002:nmae_plt_max]);
ytickformat('%1.3f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
%ax.XAxis.MinorTickValues = 1:1:num_trials;
ax.Position = [0.078 0.425 0.905 0.24];
%ax.Position = [0.10 0.44 0.88 0.23];

%--------------------------------------
subplot(3,1,3);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(3,2); num_trials+2,limits_test(3,2); num_trials+2,limits_test(3,1); -1,limits_test(3,1)];
p1 = patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);

%plot the mean line - maybe
p2 = plot([-1,num_trials+2],[ssim_stats(2,2), ssim_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(3,2), limits_test(3,2)], '-b', 'LineWidth', 1);

% plot the test points
p3 = scatter([1:num_trials], ssim(:,2),15,'filled','k');
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
% xlabel_name = {''};
% for idx=1:num_trials
%     xlabel_name{idx+1,1} = num2str(idx);
% end
xlim([0 num_trials+1]);
xticks([1:1:num_trials]);
xtickangle(90);
%xticklabels(xlabel_name(0:2:num_trials));
%set(gca,'XMinorTick','on', 'XMinorGrid','on');

xlabel(strcat('Training Event Number'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
p = 100;
ssim_plt_max = max(ceil(limits_test(3,2)*p)/p, ceil(ssim_stats(3,2)*p)/p);
ssim_plt_min = min(floor(limits_test(3,1)*p)/p, ceil(ssim_stats(1,2)*p)/p);
% ssim_plt_max = 0.925;
% ssim_plt_min = 0.895;
ylim([ssim_plt_min ssim_plt_max]);
yticks([ssim_plt_min:0.01:ssim_plt_max]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
ytickformat('%1.3f')
ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

legend([p1, p2, p3], strcat(32,'{\pm}',num2str(sig_m),'{\sigma} Metric Bounds'), 'Metric Mean', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
%ax.XAxis.MinorTickValues = 1:1:num_trials;
ax.Position = [0.078 0.145 0.905 0.24];
%ax.Position = [0.1 0.17 0.88 0.23];

if(save_results)
    print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_trial_results.png')));
    %print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_rw4_trial_results.png')));
end

plot_num = plot_num + 1;  

%% Here we will create an are plot bounded by the standard deviation
% test_label = {' ','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x'};
% test_data = [0.019033, 0.064118, 0.903725;... b
%              0.019078, 0.063447, 0.906767;... c
%              0.020007, 0.067290, 0.900115;... d
%              0.018083, 0.059135, 0.904168;... e
%              0.020541, 0.069911, 0.895514;... f
%              0.021418, 0.073557, 0.888228;... g 
%              0.019349, 0.065146, 0.900298;... h
%              0.018665, 0.064446, 0.903069;... i
%              0.022934, 0.072205, 0.879914;... j
%              0.023551, 0.074528, 0.878680;... k
%              0.022555, 0.074145, 0.884138;... l
%              0.022009, 0.073155, 0.886147;... m
%              0.019353, 0.064927, 0.896862;... n
%              0.020654, 0.068857, 0.890060;... o
%              0.022690, 0.080658, 0.884842;... p
%              0.020899, 0.067424, 0.891878;... q
%              0.020081, 0.069205, 0.893784;... r
%              0.019747, 0.069915, 0.897603;... s
%              0.019226, 0.064573, 0.900414;... t
%              0.018311, 0.061787, 0.902097;... u
%              0.019942, 0.065348, 0.895294;... v
%              0.022365, 0.072076, 0.881252;... w
%              0.018446, 0.058445, 0.902611;... x
%              ];

%% run any test data

if(run_test_data)
    
    test_listing = dir(strcat(test_data_path,filesep,'*.txt'));
    num_test_trials = numel(test_listing);

    test_label = cell(num_test_trials,1);
    test_min = zeros(num_test_trials,3);
    test_mean = zeros(num_test_trials,3);
    test_max = zeros(num_test_trials,3);
    test_std = zeros(num_test_trials,3);
    
    for idx=1:num_test_trials

        [test_label{idx,1}, ~, te_nmae, te_nrmse, te_ssim, te_nmae_stats, te_nrmse_stats, te_ssim_stats] = get_dfd_results(fullfile(test_listing(idx).folder, test_listing(idx).name), dfd_dnn_index);     

        test_min(idx,:) = [te_nmae_stats(1,2), te_nrmse_stats(1,2), te_ssim_stats(1,2)];
        test_mean(idx,:) = [te_nmae_stats(2,2), te_nrmse_stats(2,2), te_ssim_stats(2,2)];
        test_max(idx,:) = [te_nmae_stats(3,2), te_nrmse_stats(3,2), te_ssim_stats(3,2)];
        test_std(idx,:) = [te_nmae_stats(4,2), te_nrmse_stats(4,2), te_ssim_stats(4,2)];

        fprintf('\n\nName:, %s\n', test_label{idx,1});
        fprintf(' , NMAE, NRMSE, SSIM\n');
        fprintf('Min:, %2.4f, %2.4f, %2.4f\n', test_min(idx,:));
        fprintf('Mean:, %2.4f, %2.4f, %2.4f\n', test_mean(idx,:));
        fprintf('Max:, %2.4f, %2.4f, %2.4f\n', test_max(idx,:));
        fprintf('Std:, %2.4f, %2.4f, %2.4f\n', test_std(idx,:));
    end
    
    % get the overall min and max values
    nmae_min = min(min(test_min(:,1)), nmae_stats(1,2));
    nrmse_min = min(min(test_min(:,2)), nrmse_stats(1,2));
    ssim_min = min(min(test_min(:,3)), ssim_stats(1,2));
    
    nmae_max = max(max(test_max(:,1)), nmae_stats(3,2));
    nrmse_max = max(max(test_max(:,2)), nrmse_stats(3,2));
    ssim_max = max(max(test_max(:,3)), ssim_stats(3,2));    
    

%% plot the results

    small = false;
    
    min_max_dot = 13;
    error_bar_dot = 14;    

    figure(plot_num)
    if(small)
        set(gcf,'position',([50,50,800,600]),'color','w');
    else
        set(gcf,'position',([50,50,1200,700]),'color','w');
    end
    
    % Plot the NRMSE values
    subplot(3,1,1);
    hold on
    box on
    grid on

    %plot the area
    a1 = [-1,limits_test(2,2); num_trials+2,limits_test(2,2); num_trials+2,limits_test(2,1); -1,limits_test(2,1)];
    patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
    
    %plot the mean line - maybe
    plot([-1,num_trials+1],[nrmse_stats(2,2), nrmse_stats(2,2)], '--b', 'LineWidth', 1);
    plot([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], '-b', 'LineWidth', 1);
    
    % plot the test points
    errorbar([1:num_test_trials], test_mean(:,2), test_sig_m*test_std(:,2),'.',...
            'MarkerSize', error_bar_dot, 'MarkerEdgeColor','red','MarkerFaceColor','red','Color','k','LineWidth',1);    
    scatter([1:num_test_trials], test_max(:,2), min_max_dot, 'MarkerEdgeColor','green', 'MarkerFaceColor','green','LineWidth',1);
    scatter([1:num_test_trials], test_min(:,2), min_max_dot, 'MarkerEdgeColor','green', 'MarkerFaceColor','green','LineWidth',1); 
    
    set(gca, 'fontweight', 'bold', 'FontSize', 13);

    % X-Axis
    xlim([0 num_test_trials+1]);
    xticks([0:num_test_trials+1]);
    xticklabels([]);

    % Y-Axis
    plt_step = 100;
    nrmse_plt_min = floor(nrmse_min*plt_step)/plt_step;
    nrmse_plt_max = ceil(nrmse_max*plt_step)/plt_step;
%     nrmse_plt_max = 0.09;
%     nrmse_plt_min = 0.01;
    ylim([nrmse_plt_min nrmse_plt_max]);
    yticks([nrmse_plt_min:0.01:nrmse_plt_max]);
    set(gca,'YMinorTick','on', 'YMinorGrid','on');
    ytickformat('%1.3f');
    ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

    title('DFD-Net Training w/Noise Test Results','fontweight','bold','FontSize', 16)

    ax = gca;
    ax.YAxis.MinorTickValues = nrmse_plt_min:0.005:nrmse_plt_max;
    if(small)
        ax.Position = [0.1 0.725 0.88 0.21];
    else
        ax.Position = [0.07 0.715 0.92 0.22];
    end
    %-------------------------------------------------------------------------------
    subplot(3,1,2);
    hold on
    box on
    grid on

    %plot the area
    a1 = [-1,limits_test(1,2); num_trials+2,limits_test(1,2); num_trials+2,limits_test(1,1); -1,limits_test(1,1)];
    patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
    
    %plot the mean line - maybe
    plot([-1,num_trials+1],[nmae_stats(2,2), nmae_stats(2,2)], '--b', 'LineWidth', 1);
    plot([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], '-b', 'LineWidth', 1);
    
    % plot the test points
    errorbar([1:num_test_trials], test_mean(:,1), test_sig_m*test_std(:,1),'.',...
            'MarkerSize', error_bar_dot, 'MarkerEdgeColor','red','MarkerFaceColor','red','Color','k','LineWidth',1);
    scatter([1:num_test_trials], test_max(:,1), min_max_dot, 'MarkerEdgeColor','green', 'MarkerFaceColor','green','LineWidth',1);
    scatter([1:num_test_trials], test_min(:,1), min_max_dot, 'MarkerEdgeColor','green', 'MarkerFaceColor','green','LineWidth',1);
    
    set(gca, 'fontweight', 'bold', 'FontSize', 13);

    % X-Axis
    xlim([0 num_test_trials+1]);
    xticks([0:num_test_trials+1]);
    xticklabels([]);

    % Y-Axis
    plt_step = 200;
    nmae_plt_min = floor(nmae_min*plt_step)/plt_step;
    nmae_plt_max = ceil(nmae_max*plt_step)/plt_step;
    ylim([nmae_plt_min nmae_plt_max]);
    yticks([nmae_plt_min:0.010:nmae_plt_max]);
    set(gca,'YMinorTick','on', 'YMinorGrid','on');
    ytickformat('%1.3f');
    ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

    ax = gca;
    ax.YAxis.MinorTickValues = nmae_plt_min:0.002:nmae_plt_max;
    if(small)
        ax.Position = [0.1 0.48 0.88 0.21];
    else
        ax.Position = [0.07 0.465 0.92 0.22];
    end
    
    %--------------------------------------
    subplot(3,1,3);
    hold on
    box on
    grid on

    %plot the area
    a1 = [-1,limits_test(3,2); num_trials+2,limits_test(3,2); num_trials+2,limits_test(3,1); -1,limits_test(3,1)];
    p1 = patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
    
    %plot the mean line - maybe
    p2 = plot([-1,num_trials+1],[ssim_stats(2,2), ssim_stats(2,2)], '--b', 'LineWidth', 1);
    plot([-1,num_trials+2],[limits_test(3,2), limits_test(3,2)], '-b', 'LineWidth', 1);
    
    % plot the test points
    p3 = errorbar([1:num_test_trials], test_mean(:,3), test_sig_m*test_std(:,3),'.',...
        'MarkerSize', error_bar_dot, 'MarkerEdgeColor','red','MarkerFaceColor','red','Color','k','LineWidth',1);
    
    scatter([1:num_test_trials], test_max(:,3), min_max_dot, 'MarkerEdgeColor','green', 'MarkerFaceColor','green','LineWidth',1);
    scatter([1:num_test_trials], test_min(:,3), min_max_dot, 'MarkerEdgeColor','green', 'MarkerFaceColor','green','LineWidth',1);
    
    set(gca, 'fontweight', 'bold', 'FontSize', 13);

    % X-Axis
    set(gca,'TickLabelInterpreter','none')
    xlim([0 num_test_trials+1]);
    xticks([1:num_test_trials+1]);
    xticklabels(test_label);
    xlabel(strcat('Noise Type'), 'fontweight', 'bold', 'FontSize', 13);
    xtickangle(45);

    % Y-Axis
    plt_step = 20;
    ssim_plt_min = floor(ssim_min*plt_step)/plt_step;
    ssim_plt_max = ceil(ssim_max*plt_step)/plt_step;    
%     ssim_plt_max = 0.925;
%     ssim_plt_min = 0.865;
    ylim([ssim_plt_min ssim_plt_max]);
    yticks([ssim_plt_min:0.05:ssim_plt_max]);
    %ylim([0.85 1]);
    %yticks([0:0.05:1]);
    set(gca,'YMinorTick','on', 'YMinorGrid','on');
    ytickformat('%1.3f');
    ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

    legend([p1, p2, p3], strcat('{\pm}',32,num2str(sig_m),'{\sigma} Metric Bounds'), 'Metric Mean', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')
    ax = gca;
    ax.YAxis.MinorTickValues = ssim_plt_min:0.010:ssim_plt_max;
    if(small)
        ax.Position = [0.1 0.235 0.88 0.21];
    else
        ax.Position = [0.07 0.21 0.92 0.22];
    end
    
    if(save_results)
        print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_noise_results.png')));
    end
    
    plot_num = plot_num + 1;  
 
end

%% combined table to show both training and testing results
fprintf('\n\n');

caption = {'DFD-Net Average Training Performance Results', 'DFD-Net Average Test Performance Results'};
tbl_name = {'tbl:dnn_trials_1', 'tbl:dnn_trial_2'};
row_name = {'Minimum', 'Mean', 'Maximum', 'Std Deviation'};

write_latex_table_head('DFD-Net Average Training/Testing Performance Results', 'tbl:dnn_trials1', 'l|c|c|c|c|c|c|')
fprintf('\\cline{2-7} \n');
fprintf('\\multicolumn{1}{c|}{} & \\multicolumn{3}{c|}{\\textbf{Training}} & \\multicolumn{3}{c|}{\\textbf{Testing}} \\\\ \\cline{2-7} \n');
fprintf('\\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{\\textbf{NRMSE}} & \\multicolumn{1}{c|}{\\textbf{NMAE}} & \\multicolumn{1}{c|}{\\textbf{SSIM}} & \\multicolumn{1}{c|}{\\textbf{NRMSE}} & \\multicolumn{1}{c|}{\\textbf{NMAE}} & \\multicolumn{1}{c|}{\\textbf{SSIM}} \\\\ \\hline \n');
% print out the values
for idx=1:4
    fprintf('\\multicolumn{1}{|l|}{\\textbf{%s}} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', row_name{idx}, nrmse_stats(idx,1), nmae_stats(idx,1), ssim_stats(idx,1), nrmse_stats(idx,2), nmae_stats(idx,2), ssim_stats(idx,2));   
end
write_latex_table_end;

fprintf('\n\n');

