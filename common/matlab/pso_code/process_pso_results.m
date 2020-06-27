format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% Select PSO results files

save_path = 'D:\IUPUI\PhD\Images\PSO';

file_filter = {'*.txt','Text Files';'*.*','All Files' };
startpath = 'D:\IUPUI\DfD\dfd_dnn_pso\pso_results\';

[results_file, results_file_path] = uigetfile(file_filter, 'Select the PSO Results File(s)', startpath, 'multiselect','on');
if(results_file_path == 0)
    return;
end

save_results = false;

save_path = 'D:\IUPUI\PhD\Images\PSO';
save_path = uigetdir(save_path,'Select Save Folder');

if(save_path ~= 0)
    save_results = true;
end


commandwindow;

%%  Process through the files and read in the data

N = 20;                         % population size
gamma = 0.3;
cf = 255/224;

if(iscell(results_file))
    
    itr = numel(results_file);

    % arrange the data in column form where the PSO member is the column
    results = zeros(N, 6, itr);
    f = zeros(N,itr);
    
    for idx=1:itr
        tmp = parse_input_parameters(fullfile(results_file_path, results_file{idx}));
        
        for jdx=1:numel(tmp)
        
            % read in the iteration number
            itr_num = str2double(tmp{jdx}{1});
            
            % read in the population member
            pop_num = str2double(tmp{jdx}{2});
        
            results(pop_num,1,itr_num) = str2double(tmp{jdx}{3})*cf;
            results(pop_num,2,itr_num) = str2double(tmp{jdx}{4})*cf;
            results(pop_num,3,itr_num) = str2double(tmp{jdx}{5});
            results(pop_num,4,itr_num) = str2double(tmp{jdx}{8})*cf;
            results(pop_num,5,itr_num) = str2double(tmp{jdx}{9})*cf;
            results(pop_num,6,itr_num) = str2double(tmp{jdx}{10});
            
            f(pop_num,itr_num) = gamma*(results(pop_num,1,itr_num) + results(pop_num,2,itr_num) + (1-results(pop_num,3,itr_num))) + ...
                                (1-gamma)*(results(pop_num,4,itr_num) + results(pop_num,5,itr_num) + (1-results(pop_num,6,itr_num)));
            
        end
    end
    
else
    itr = 1;
    results = zeros(N, 6, itr);

end

itr_num = itr;

%% get some of the basic metrics

particle_min = zeros(N,6);
particle_mean = zeros(N,6);
particle_max= zeros(N,6);

for idx=1:N
    particle_min(idx,:) = min(results(idx,:,:),[],3);
    particle_mean(idx,:) = mean(results(idx,:,:),3);
    particle_max(idx,:) = max(results(idx,:,:),[],3);
end

min_metrics = zeros(6,3);
mean_metrics = zeros(6,3);
max_metrics = zeros(6,3);

for idx=1:6
    
    [tmp, index_x] = min(results(:,idx,:));
    [tmp, index_y] = min(tmp);
    min_metrics(idx,:) = [index_x(:,:,index_y),index_y, tmp];
    
%     [tmp, index_x] = mean(results(:,idx,:));
%     [tmp, index_y] = mean(tmp);
%     mean_metrics(idx,:) = [index_x(:,:,index_y),index_y, tmp];
    
    [tmp, index_x] = max(results(:,idx,:));
    [tmp, index_y] = max(tmp);
    max_metrics(idx,:) = [index_x(:,:,index_y),index_y, tmp];
end

best_result = [min_metrics(1,:); min_metrics(2,:); max_metrics(3,:); min_metrics(4,:); min_metrics(5,:); max_metrics(6,:)];


f_min = min(f,[],1);
f_max = max(f,[],1);
f_mean = mean(f,1);


%% plot the objective function stats

figure(plot_num)
set(gcf,'position',([100,50,1200,700]),'color','w')
hold on
box on
grid on

p1 = plot([1:itr], f_min, 'g', 'LineWidth', 1);
p2 = plot([1:itr], f_mean, '--k', 'LineWidth', 1);
p3 = plot([1:itr], f_max, 'b', 'LineWidth', 1);

set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([1 itr]);
xticks([1:itr]);
xlabel('PSO Iteration Number','fontweight', 'bold', 'FontSize', 13);

% Y-axis
ylim([0 2.1]);
yticks([0:0.1:2.1]);
ytickformat('%1.1f')
ylabel('Objective Function Value', 'fontweight', 'bold', 'FontSize', 13);

title('Objective Function', 'fontweight', 'bold', 'FontSize', 16);
lgd = legend([p1, p2, p3], 'Minimum', 'Mean', 'Maximum', 'location', 'northeast', 'orientation', 'vertical');

ax = gca;
ax.Position = [0.06 0.09 0.91 0.85];

%lgd.Position = [.36, .02, .32, 0.034];

if(save_results)
    %print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_pso_fx_results.png')));
end

plot_num = plot_num + 1;

%% setup the plot limits
% which one to plot: train: NMAE=1, NRMSE=2, SSIM=3, test: NMAE=4, NRMSE=5, SSIM=6
metric = 4;
metric_name = {'Training NMAE','Training NRMSE','Training SSIM', 'Test NMAE','Test NRMSE','Test SSIM'};
metric_save_name = {'nmae_tr', 'nrmse_tr', 'ssim_tr', 'nmae_te', 'nrmse_te', 'ssim_te'};

% z is the actual value [nmae_min,nmae_max;nrmse_min,nrmse_max;ssim_min,ssim_max]
z_lim = [0.00,0.60; 0.00,0.6; 0.00,1.00; ...
         0.00,0.60; 0.00,0.6; 0.00,1.00];
     
z_step = [0.05,0.05,0.2, 0.05,0.05,0.2];

z_ticks = [z_lim(metric,1):z_step(metric):z_lim(metric,2)];


%% 3-D view of the scatter plot

figure(plot_num)
set(gcf,'position',([100,50,1200,700]),'color','w')
hold on
box on
grid on

% Y = 1:1:N;
% Y = repmat(Y,1,itr_num);
% 
% X = 1:1:itr_num;
% X = repmat(X,1,N);
% 
% Z = permute(results(:,metric,:),[3 2 1]);
% Z = Z(:)';

S = 15*ones(1,itr_num);
cm = colormap(jet(ceil(z_lim(metric,2)*100))); 

for idx=1:N
    
    X = idx*ones(1,itr_num);    % this will be the population member
    Y = 1:itr_num;    % this will be the iteration number
    Z = results(idx,metric,:);    % this will be the data
                                        % Define Colormap
    c2 = ceil(Z(:)*100);
    
    s1 = scatter3(X, Y, Z, 15, 'filled', 'b');
    %s1 = scatter3(X, Y, Z(:), S, cm(c2,:), 'filled');
end

sb = scatter3(best_result(metric,1), best_result(metric,2), best_result(metric,3), 20, 'filled', 'r');

%s1 = scatter3(X, Y, Z, 15, 'filled', 'b');
%sb = scatter3(best_result(metric,2), best_result(metric,1), best_result(metric,3), 18, 'filled', 'r');

set(gca,'fontweight','bold');

% Y-Axis
ylim([-1 itr_num]);
yticks([0:1:itr_num-1]);
%ytickformat('%1.2f')
ylabel(strcat('Iteration Number',13), 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
%set(gca,'TickLabelInterpreter','none')
xlim([0 N+1]);
xticks([1:1:N]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xlabel(strcat('Particle Number',13), 'fontweight', 'bold', 'FontSize', 13);
%xticklabels(test_label);

% Z-axis
zlim(z_lim(metric,:))
zticks(z_ticks);
zlabel(metric_name{metric}, 'fontweight', 'bold', 'FontSize', 13);
ztickformat('%1.2f')
view(60,45);

%title('MNIST Multi-Training Event - Test Results Distribution','fontweight','bold','FontSize', 15)

lgd = legend([s1, sb], 'Particle Results', 'Global Best Particle', 'location', 'southoutside', 'orientation', 'horizontal');

ax = gca;
ax.Position = [0.08 0.15 0.88 0.82];

lgd.Position = [.36, .02, .32, 0.034];

if(save_results)
    %print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_pso_',metric_save_name{metric},'_results.png')));
end

plot_num = plot_num + 1;


%% surface version of the plot over iteration and population member

if(false)
% build the surface
X = 1:N;
Y = 1:itr_num;

Z = zeros(itr_num,N);

for idx=1:N
   Z(:,idx) = results(idx,metric,:);
end

figure(plot_num)
set(gcf,'position',([100,100,1000,600]),'color','w')
hold on
box on
grid on

s1 = surf(X,Y,Z);
colormap(jet(1000));
sb = scatter3(best_result(metric,1), best_result(metric,2), best_result(metric,3), 18, 'filled', 'k');

set(gca,'fontweight','bold');

% Y-Axis
ylim([-1 itr_num]);
yticks([0:1:itr_num-1]);
%ytickformat('%1.2f')
ylabel('Iteration Number', 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
%set(gca,'TickLabelInterpreter','none')
xlim([0 N+1]);
xticks([1:1:N]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xlabel('Particle Number', 'fontweight', 'bold', 'FontSize', 13);
%xticklabels(test_label);

% Z-axis
zlim(z_lim(metric,:))
zticks(z_ticks);
zlabel(metric_name{metric}, 'fontweight', 'bold', 'FontSize', 13);

view(60,45);

plot_num = plot_num + 1;

end
%% 3-D scatter plot with both train and test

if(false)

m2 = mod(metric,3)+1;
    
z_label_name = {'NMAE','NRMSE','SSIM'};

figure(plot_num)
set(gcf,'position',([100,100,1000,600]),'color','w')
hold on
box on
grid on

for idx=1:N

    % this will be the population member
    X = idx*ones(itr_num,1);

    % this will be the iteration number
    Y = 1:itr_num;

    % this will be the data
    Z1 = results(idx,m2,:);
    Z2 = results(idx,m2+3,:);

    s1 = scatter3(X, Y, Z1, 15, 'filled', 'b');   
    s2 = scatter3(X, Y, Z2, 15, 'filled', 'g');

end

set(gca,'fontweight','bold');

% Y-Axis
ylim([-1 itr_num]);
yticks([0:1:itr_num-1]);
%ylim([0.85 1]);
%yticks([0:0.05:1]);
%ytickformat('%1.2f')
ylabel('Iteration Number', 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
%set(gca,'TickLabelInterpreter','none')
xlim([0 N+1]);
xticks([1:1:N]);
%set(gca,'XMinorTick','on', 'XMinorGrid','on');
xlabel('Particle Number', 'fontweight', 'bold', 'FontSize', 13);
%xticklabels(test_label);

% Z-axis
zlim(z_lim(metric,:))
zticks(z_ticks);
zlabel(z_label_name{m2}, 'fontweight', 'bold', 'FontSize', 13);

view(60,45);

%title('MNIST Multi-Training Event - Test Results Distribution','fontweight','bold','FontSize', 15)

plot_num = plot_num + 1;

end

%% plot the particle results on the training distribution graph
num_trials = N;
sig_m = 2;

nrmse_stats = [0.0428,0.0552; 0.0473,0.0631; 0.0497,0.0699; 0.0012,0.0039];
nmae_stats = [0.0104,0.0163; 0.0128,0.0183; 0.0137,0.0200; 0.0006,0.0009];
ssim_stats = [0.9354,0.9000; 0.9391,0.9074; 0.9469,0.9180; 0.0022,0.0044];

% get the lower limits
limits_test(1,1) = nmae_stats(2,2) - sig_m*nmae_stats(4,2);
limits_test(2,1) = nrmse_stats(2,2) - sig_m*nrmse_stats(4,2);
limits_test(3,1) = ssim_stats(2,2) - sig_m*ssim_stats(4,2);

% get the upper limits
limits_test(1,2) = nmae_stats(2,2) + sig_m*nmae_stats(4,2);
limits_test(2,2) = nrmse_stats(2,2) + sig_m*nrmse_stats(4,2);
limits_test(3,2) = ssim_stats(2,2) + sig_m*ssim_stats(4,2);

[g_best_nrmse, g_best_nrmse_index] = min(particle_min(:,5));
[g_best_nmae, g_best_nmae_index] = min(particle_min(:,4));
[g_best_ssim, g_best_ssim_index] = max(particle_max(:,6));

figure(plot_num)
set(gcf,'position',([50,50,1200,700]),'color','w')
%set(gcf,'position',([50,50,1000,600]),'color','w')

% Plot the NRMSE values
subplot(3,1,1);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(2,2); num_trials+2,limits_test(2,2); num_trials+2,limits_test(2,1); -1,limits_test(2,1)];
%a2 = [-1,limits_test(5,2); num_scene+2,limits_test(5,2); num_scene+2,limits_test(5,1); -1,limits_test(5,1)];

patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%patch('Faces',[1 2 3 4], 'Vertices',a2,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%area([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], 'BaseValue', limits_test(2,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth', 1)

%plot the mean line - maybe
plot([-1,num_trials+2],[nrmse_stats(2,2), nrmse_stats(2,2)], '--b', 'LineWidth', 1);
%plot([-1,num_trials+2],[limits_test(2,2), limits_test(2,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], particle_min(:,5), 17, 'filled', 'k');
scatter(g_best_nrmse_index, g_best_nrmse, 20, 'filled', 'g');

set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
set(gca,'XMinorTick','on', 'XMinorGrid','on');
xticklabels([]);

% Y-Axis
%nrmse_plt_max = ceil(limits_test(2,2)*80)/80;
%nrmse_plt_min = floor(limits_test(2,1)*80)/80;
nrmse_plt_max = 0.08;
nrmse_plt_min = 0.04;
ylim([nrmse_plt_min nrmse_plt_max]);
yticks([nrmse_plt_min:0.01:nrmse_plt_max]);
ytickformat('%1.3f')
ylabel('NRMSE', 'fontweight', 'bold', 'FontSize', 13);

title('DfD-Net PSO Particle Best Results','fontweight','bold','FontSize', 16)

ax = gca;
ax.XAxis.MinorTickValues = 1:1:num_trials;
%ax.Position = [0.078 0.700 0.905 0.24];
ax.Position = [0.07 0.710 0.92 0.24];

%--------------------------------------
subplot(3,1,2);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(1,2); num_trials+2,limits_test(1,2); num_trials+2,limits_test(1,1); -1,limits_test(1,1)];
%a2 = [-1,limits_test(4,2); num_scene+2,limits_test(4,2); num_scene+2,limits_test(4,1); -1,limits_test(4,1)];

patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%patch('Faces',[1 2 3 4], 'Vertices',a2,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%area([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], 'BaseValue', limits_test(1,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1)

%plot the mean line - maybe
plot([-1,num_trials+2],[nmae_stats(2,2), nmae_stats(2,2)], '--b', 'LineWidth', 1);
%plot([-1,num_trials+2],[limits_test(1,2), limits_test(1,2)], '-b', 'LineWidth', 1);

% plot the test points
scatter([1:num_trials], particle_min(:,4), 17, 'filled', 'k');
scatter(g_best_nmae_index, g_best_nmae, 20, 'filled', 'g');

set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
xlim([0 num_trials+1]);
xticks([0:num_trials+1]);
set(gca,'XMinorTick','on', 'XMinorGrid','on');
xticklabels([]);

% Y-Axis
%nmae_plt_max = ceil(limits_test(1,2)*200)/200;
%nmae_plt_min = floor(limits_test(1,1)*200)/200;
nmae_plt_max = 0.022;
nmae_plt_min = 0.014;
ylim([nmae_plt_min nmae_plt_max]);
yticks([nmae_plt_min:0.002:nmae_plt_max]);
ytickformat('%1.3f')
ylabel('NMAE', 'fontweight', 'bold', 'FontSize', 13);

ax = gca;
ax.XAxis.MinorTickValues = 1:1:num_trials;
%ax.Position = [0.078 0.425 0.905 0.24];
ax.Position = [0.07 0.425 0.92 0.24];

%--------------------------------------
subplot(3,1,3);
hold on
box on
grid on

%plot the area
a1 = [-1,limits_test(3,2); num_trials+2,limits_test(3,2); num_trials+2,limits_test(3,1); -1,limits_test(3,1)];
%a2 = [-1,limits_test(6,2); num_scene+2,limits_test(6,2); num_scene+2,limits_test(6,1); -1,limits_test(6,1)];

p1=patch('Faces',[1 2 3 4], 'Vertices',a1,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%patch('Faces',[1 2 3 4], 'Vertices',a2,'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
%p1 = area([-1,num_trials+2],[limits_test(3,2), limits_test(3,2)], 'BaseValue', limits_test(3,1), 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',1, 'EdgeColor','b', 'LineWidth',1);

%plot the mean line - maybe
p2 = plot([-1,num_trials+2],[ssim_stats(2,2), ssim_stats(2,2)], '--b', 'LineWidth', 1);
%plot([-1,num_trials+2],[limits_test(3,2), limits_test(3,2)], '-b', 'LineWidth', 1);

% plot the test points
p3 = scatter([1:num_trials], particle_max(:,6),17,'filled','k');
p4 = scatter(g_best_ssim_index, g_best_ssim, 20, 'filled', 'g');
set(gca, 'fontweight', 'bold', 'FontSize', 13);

% X-Axis
% xlabel_name = {''};
% for idx=1:num_trials
%     xlabel_name{idx+1,1} = num2str(idx);
% end
xlim([0 num_trials+1]);
xticks([1:1:num_trials]);
%xticklabels(xlabel_name(0:2:num_trials));
%set(gca,'XMinorTick','on', 'XMinorGrid','on');

xlabel(strcat('Particle Number'), 'fontweight', 'bold', 'FontSize', 13);

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

legend([p1, p2, p3, p4], strcat(32,'{\pm}',num2str(sig_m),'{\sigma} Metric Bounds'), 'Metric Mean', 'P-Best Results', 'G-Best Results', 'location', 'southoutside', 'orientation', 'horizontal')
ax = gca;
ax.XAxis.MinorTickValues = 1:1:num_trials;
%ax.Position = [0.078 0.145 0.905 0.24];
ax.Position = [0.07 0.14 0.92 0.24];

if(save_results)
    % print(plot_num, '-dpng', fullfile(save_path,strcat('dfd_pso_best_results.png')));
end

plot_num = plot_num + 1;  


%% The end
fprintf('Complete!\n');


