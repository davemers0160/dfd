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

%% Get the data from the file

data_params = parse_input_parameters(fullfile(data_path, data_file));

num_trials = size(data_params,1);

nmae = zeros(num_trials,2);
nrmse = zeros(num_trials,2);
ssim = zeros(num_trials,2);

nmae = zeros(num_trials,1);
nrmse = zeros(num_trials,1);
ssim = zeros(num_trials,1);

for idx=1:num_trials
    nmae(idx,1) = str2double(data_params{idx}{1});
    nrmse(idx,1) = str2double(data_params{idx}{2});
    ssim(idx,1) = str2double(data_params{idx}{3});
    nmae(idx,2) = str2double(data_params{idx}{4});
    nrmse(idx,2) = str2double(data_params{idx}{5});
    ssim(idx,2) = str2double(data_params{idx}{6}); 
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

nrmse_plt_max = ceil(max(nrmse(:))*25)/25;
nmae_plt_max = ceil(max(nmae(:))*34)/34;

lw = 1.5;
sig_m = 2;

%% plot the training results 
if(false)
    
    title_str = {'DFD-Net Average Training Performance Results', 'DFD-Net Average Test Performance Results'};
    plot_save_name = {'_v14a_train', '_v14a_test'};

    for idx=1:2

        fig = figure(plot_num);
        set(gcf,'position',([50,50,1000,700]),'color','w')

        % Plot the NRMSE values
        subplot(3,1,1);
        hold on
        box on
        grid on

        plot([1:num_trials], nrmse(:,idx) ,'-b');
        set(gca, 'fontweight', 'bold', 'FontSize', 13);

        % X-Axis
        xlim([1 num_trials]);
        xticks([1:num_trials]);
        xticklabels([]);

        % Y-Axis
        ylim([0 nrmse_plt_max]);
        yticks([0:0.02:nrmse_plt_max]);
        ytickformat('%1.2f')
        ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

        title(strcat(title_str{idx}),'fontweight','bold','FontSize', 16)

        ax = gca;
        ax.Position = [0.07 0.675 0.915 0.245];

        %--------------------------------------
        % Plot the NMAE values
        subplot(3,1,2);
        hold on
        box on
        grid on

        plot([1:num_trials], nmae(:,idx) ,'-b');
        set(gca,'fontweight','bold','FontSize', 13);

        % X-Axis
        xlim([1 num_trials]);
        xticks([1:num_trials]);
        xticklabels([]);

        % Y-Axis
        ylim([0 nmae_plt_max]);
        yticks([0:0.01:nmae_plt_max]);
        ytickformat('%1.2f');
        ylabel('NMAE','fontweight','bold','FontSize', 13);

        %title(strcat('DFD-Net Average Training Performance Results'),'fontweight','bold','FontSize', 16)

        ax = gca;
        ax.Position = [0.07 0.395 0.915 0.245];

        %--------------------------------------
        % Plot the SSIM values
        subplot(3,1,3);
        hold on
        box on
        grid on

        plot([1:num_trials], ssim(:,idx) ,'-b');
        set(gca,'fontweight','bold','FontSize', 13);

        % X-Axis
        xlim([1 num_trials]);
        xticks([1:num_trials]);
        %xtickangle(90);
        %xtickformat();
        xlabel(strcat('Training Trial Number'), 'fontweight', 'bold', 'FontSize', 13);

        % Y-Axis
        ylim([0 1]);
        yticks([0:0.25:1]);
        ytickformat('%1.2f');
        ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

        ax = gca;
        ax.Position = [0.07 0.11 0.915 0.245];

        print(plot_num, '-dpng', fullfile(data_path,strcat('dfd_dnn_perf',plot_save_name{idx},'.png')));

        plot_num = plot_num + 1;  
    end
    
end

%% combined training and testing plot

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
    yticks([0:0.02:nrmse_plt_max]);
    ytickformat('%1.2f')
    ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

    title('Multiple Trial Performance Comparison','fontweight','bold','FontSize', 16)

    ax = gca;
    ax.Position = [0.07 0.700 0.915 0.24];

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
    yticks([0:0.01:nmae_plt_max]);
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

    plot([1:num_trials], ssim(:,1) ,'--b', 'LineWidth', lw);
    plot([1:num_trials], ssim(:,2) ,'-b', 'LineWidth', lw);
    set(gca,'fontweight','bold','FontSize', 13);

    % X-Axis
    xlim([1 num_trials]);
    xticks([1:num_trials]);
    %xtickangle(90);
    %xtickformat();
    xlabel(strcat('Training Trial Number'), 'fontweight', 'bold', 'FontSize', 13);

    % Y-Axis
    ylim([0.85 1]);
    yticks([0:0.05:1]);
    ytickformat('%1.2f');
    ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);
    legend('Training Results', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')

    ax = gca;
    ax.Position = [0.07 0.145 0.915 0.24];

    print(plot_num, '-dpng', fullfile(data_path,strcat('dfd_dnn_combine_mutli_trial.png')));

    plot_num = plot_num + 1;  

end

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
set(gcf,'position',([50,50,1000,700]),'color','w')
% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
area([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], 'BaseValue', limits_test(2,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth', 1)

%plot the mean line - maybe
plot([-1,num_trials+2],[nrmse_stats(2,2), nrmse_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], nrmse(:,2),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
xticklabels([]);

% Y-Axis
%nrmse_plt_max = ceil(limits_test(2,2)*80)/80;
%nrmse_plt_min = floor(limits_test(2,1)*80)/80;
nrmse_plt_max = 0.08;
nrmse_plt_min = 0.05;
ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:0.01:nrmse_plt_max]);
ytickformat('%1.3f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

title('Multiple Trial Performance Comparison','fontweight','bold','FontSize', 16)

ax = gca;
ax.Position = [0.078 0.700 0.905 0.24];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
area([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], 'BaseValue', limits_test(1,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+2],[nmae_stats(2,2), nmae_stats(2,2)], '--b', 'LineWidth', 1);
plot([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], nmae(:,2),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
xticklabels([]);

% Y-Axis
%nmae_plt_max = ceil(limits_test(1,2)*200)/200;
%nmae_plt_min = floor(limits_test(1,1)*200)/200;
nmae_plt_max = 0.021;
nmae_plt_min = 0.015;
ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:0.002:nmae_plt_max]);
ytickformat('%1.3f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
ax.Position = [0.078 0.425 0.905 0.24];
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
p3 = scatter([1:num_trials], ssim(:,2),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlabel_name = {''};
for idx=1:num_trials
    xlabel_name{idx+1,1} = num2str(idx);
end
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
xticklabels(xlabel_name);
xlabel(strcat('Trial Number'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
%ssim_plt_max = ceil(limits_test(3,2)*80)/80;
%ssim_plt_min = floor(limits_test(3,1)*80)/80;
ssim_plt_max = 0.925;
ssim_plt_min = 0.895;
ylim([ssim_plt_min ssim_plt_max]);
yticks([ssim_plt_min:0.01:ssim_plt_max]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
ytickformat('%1.3f')
ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

legend([p1, p2, p3], strcat(32,'{\pm}',num2str(sig_m),'{\sigma} Trial Bounds'), 'Trial Mean', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.078 0.145 0.905 0.24];

plot_num = plot_num + 1;  

%% Here we will create an are plot bounded by the standard deviation
test_label = {' ','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x'};
test_data = [0.019033, 0.064118, 0.903725;... b
             0.019078, 0.063447, 0.906767;... c
             0.020007, 0.067290, 0.900115;... d
             0.018083, 0.059135, 0.904168;... e
             0.020541, 0.069911, 0.895514;... f
             0.021418, 0.073557, 0.888228;... g 
             0.019349, 0.065146, 0.900298;... h
             0.018665, 0.064446, 0.903069;... i
             0.022934, 0.072205, 0.879914;... j
             0.023551, 0.074528, 0.878680;... k
             0.022555, 0.074145, 0.884138;... l
             0.022009, 0.073155, 0.886147;... m
             0.019353, 0.064927, 0.896862;... n
             0.020654, 0.068857, 0.890060;... o
             0.022690, 0.080658, 0.884842;... p
             0.020899, 0.067424, 0.891878;... q
             0.020081, 0.069205, 0.893784;... r
             0.019747, 0.069915, 0.897603;... s
             0.019226, 0.064573, 0.900414;... t
             0.018311, 0.061787, 0.902097;... u
             0.019942, 0.065348, 0.895294;... v
             0.022365, 0.072076, 0.881252;... w
             0.018446, 0.058445, 0.902611;... x
             ];

sig_m = 2;

% get the lower limits
limits_test(1,1) = nmae_stats(2,2) - sig_m*nmae_stats(4,2);
limits_test(2,1) = nrmse_stats(2,2) - sig_m*nrmse_stats(4,2);
limits_test(3,1) = ssim_stats(2,2) - sig_m*ssim_stats(4,2);

% get the upper limits
limits_test(1,2) = nmae_stats(2,2) + sig_m*nmae_stats(4,2);
limits_test(2,2) = nrmse_stats(2,2) + sig_m*nrmse_stats(4,2);
limits_test(3,2) = ssim_stats(2,2) + sig_m*ssim_stats(4,2);


figure(plot_num)
set(gcf,'position',([50,50,1000,700]),'color','w')
% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
area([-1,num_trials+1],[limits_test(2,2), limits_test(2,2)], 'BaseValue', limits_test(2,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','k', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+1],[nrmse_stats(2,2), nrmse_stats(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:size(test_data,1)],test_data(:,2),20,'filled','r')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials]);
xticks([0:num_trials]);
xticklabels([]);

% Y-Axis
nrmse_plt_max = ceil(limits_test(2,2)*35)/35;
nrmse_plt_min = floor(limits_test(2,1)*80)/80;
ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:0.01:nrmse_plt_max]);
ytickformat('%1.3f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

title('Multiple Trial Performance Comparison','fontweight','bold','FontSize', 16)

ax = gca;
ax.Position = [0.078 0.700 0.905 0.24];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
area([-1,num_trials+1],[limits_test(1,2), limits_test(1,2)], 'BaseValue', limits_test(1,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','k', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+1],[nmae_stats(2,2), nmae_stats(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:size(test_data,1)],test_data(:,1),20,'filled','r')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials]);
xticks([0:num_trials]);
xticklabels([]);

% Y-Axis
nmae_plt_max = ceil(limits_test(1,2)*200)/200;
nmae_plt_min = floor(limits_test(1,1)*200)/200;
ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:0.002:nmae_plt_max]);
ytickformat('%1.3f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
ax.Position = [0.078 0.425 0.905 0.24];
%--------------------------------------
subplot(3,1,3);
hold on
box on
grid on

%plot the area
area([-1,num_trials+1],[limits_test(3,2), limits_test(3,2)], 'BaseValue', limits_test(3,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','k', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+1],[ssim_stats(2,2), ssim_stats(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:size(test_data,1)],test_data(:,3),20,'filled','r')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials]);
xticks([0:num_trials]);
xticklabels(test_label);
xlabel(strcat('Trial Number'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
%ssim_plt_max = ceil(limits_test(3,2)*80)/80;
%ssim_plt_min = floor(limits_test(3,1)*80)/80;
ssim_plt_max = 0.925;
ssim_plt_min = 0.875;
ylim([ssim_plt_min ssim_plt_max]);
yticks([ssim_plt_min:0.01:ssim_plt_max]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
ytickformat('%1.3f')
ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

legend(strcat(num2str(sig_m),'{\sigma} Trial Bounds'), 'Trial Mean', 'Test Results', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.078 0.145 0.905 0.24];

plot_num = plot_num + 1;  

%% Plot the path size data

test_label = {' ','12x12','16x16','20x20','24x24','28x28','32x32','36x36','40x40','44x44','48x48','52x52','56x56','60x60','64x64','68x68'};

test_data = [...0.3251,0.4081,0.2366,0.2887,0.3803,0.2269;...
             0.0348,0.0784,0.8556,0.0358,0.0881,0.8228;...
             0.0124,0.0428,0.9358,0.0216,0.0649,0.8887;...
             0.0185,0.0559,0.9113,0.0227,0.0665,0.8863;...
             0.0120,0.0462,0.9405,0.0191,0.0647,0.9036;...
             0.0123,0.0458,0.9399,0.0189,0.0656,0.9051;...
             0.0133,0.0474,0.9391,0.0173,0.0600,0.9132;...
             0.0147,0.0505,0.9349,0.0184,0.0609,0.9078;...
             0.0137,0.0496,0.9351,0.0196,0.0691,0.8995;...
             0.0141,0.0502,0.9344,0.0202,0.0724,0.8960;...
             0.0156,0.0526,0.9291,0.0225,0.0771,0.8871;...
             0.0176,0.0566,0.9210,0.0238,0.0773,0.8811;...
             0.0200,0.0615,0.9130,0.0266,0.0840,0.8709;...
             0.0207,0.0606,0.9084,0.0250,0.0758,0.8768;...
             0.0215,0.0623,0.9082,0.0263,0.0779,0.8753;...
             0.0267,0.0707,0.8790,0.0265,0.0685,0.8677;...
             ];
td_size = size(test_data,1);
sig_m = 2;

% get the lower limits
limits_test(1,1) = nmae_stats(2,2) - sig_m*nmae_stats(4,2);
limits_test(2,1) = nrmse_stats(2,2) - sig_m*nrmse_stats(4,2);
limits_test(3,1) = ssim_stats(2,2) - sig_m*ssim_stats(4,2);

% get the upper limits
limits_test(1,2) = nmae_stats(2,2) + sig_m*nmae_stats(4,2);
limits_test(2,2) = nrmse_stats(2,2) + sig_m*nrmse_stats(4,2);
limits_test(3,2) = ssim_stats(2,2) + sig_m*ssim_stats(4,2);


figure(plot_num)
set(gcf,'position',([50,50,1000,700]),'color','w')
% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
area([-1,num_trials+1],[limits_test(2,2), limits_test(2,2)], 'BaseValue', limits_test(2,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','k', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+1],[nrmse_stats(2,2), nrmse_stats(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:size(test_data,1)],test_data(:,5),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 td_size+1]);
xticks([0:td_size+1]);
xticklabels([]);

% Y-Axis
nrmse_plt_max = 0.09;   %ceil(limits_test(2,2)*35)/35;
nrmse_plt_min = 0.05;   %floor(limits_test(2,1)*80)/80;
ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:0.01:nrmse_plt_max]);
ytickformat('%1.3f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

title('Multiple Trial Performance Comparison','fontweight','bold','FontSize', 16)

ax = gca;
ax.Position = [0.078 0.705 0.905 0.235];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
area([-1,num_trials+1],[limits_test(1,2), limits_test(1,2)], 'BaseValue', limits_test(1,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','k', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+1],[nmae_stats(2,2), nmae_stats(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:size(test_data,1)],test_data(:,4),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 td_size+1]);
xticks([0:td_size+1]);
xticklabels([]);

% Y-Axis
nmae_plt_max = 0.036;   %ceil(limits_test(1,2)*200)/200;
nmae_plt_min = 0.015;   %floor(limits_test(1,1)*200)/200;
ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:0.005:nmae_plt_max]);
ytickformat('%1.3f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
ax.Position = [0.078 0.435 0.905 0.235];
%--------------------------------------
subplot(3,1,3);
hold on
box on
grid on

%plot the area
area([-1,num_trials+1],[limits_test(3,2), limits_test(3,2)], 'BaseValue', limits_test(3,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','k', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+1],[ssim_stats(2,2), ssim_stats(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:size(test_data,1)],test_data(:,6),15,'filled','k')
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 td_size+1]);
xticks([0:td_size+1]);
xticklabels(test_label);
xtickangle(90);
xlabel(strcat('Trial Number'), 'fontweight', 'bold', 'FontSize', 13);

% Y-Axis
%ssim_plt_max = ceil(limits_test(3,2)*80)/80;
%ssim_plt_min = floor(limits_test(3,1)*80)/80;
ssim_plt_max = 0.92;
ssim_plt_min = 0.82;
ylim([ssim_plt_min ssim_plt_max]);
yticks([ssim_plt_min:0.02:ssim_plt_max]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
ytickformat('%1.3f')
ylabel('SSIM', 'fontweight', 'bold', 'FontSize', 13);

legend(strcat(num2str(sig_m),'{\sigma} Trial Bounds'), 'Trial Mean', 'Patch Size', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.Position = [0.078 0.16 0.905 0.235];

plot_num = plot_num + 1;  

%% print out the latex tables

fprintf('\n\n');

caption = {'DFD-Net Average Training Performance Results', 'DFD-Net Average Test Performance Results'};
tbl_name = {'tbl:dnn_trials_1', 'tbl:dnn_trial_2'};
row_name = {'minimum', 'mean', 'maximum', 'std dev'};

for jdx=1:2

    write_latex_table_head(caption{jdx}, tbl_name{jdx}, 'l|c|c|c|')
    fprintf('\\cline{2-4} \n');
    fprintf('\\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{\\textbf{NRMSE}} & \\multicolumn{1}{c|}{\\textbf{NMAE}} & \\multicolumn{1}{c|}{\\textbf{SSIM}} \\\\ \\hline \n');
    % print out the values
    for idx=1:4
        fprintf('\\multicolumn{1}{|l|}{%s} & %2.4f & %2.4f & %2.4f \\\\ \\hline \n', row_name{idx}, nrmse_stats(idx,jdx), nmae_stats(idx,jdx), ssim_stats(idx,jdx));   
    end
    write_latex_table_end;

    fprintf('\n\n');

end

%% combined table t show both training and testing results
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

