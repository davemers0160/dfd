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

%% get the log file to parse
if(ispc)
    startpath = 'D:\IUPUI\PhD\Results';
else
    startpath = '/home/owner/DfD';
end

[log_file, log_path] = uigetfile(file_filter, 'Select DFD-Net Log File', startpath);
if(log_path == 0)
    return;
end

%% Start reading the log file
[log_version, data, final_traning, final_test] = read_dnn_logfile(fullfile(log_path,log_file));


%% 
names = fieldnames(data);

[lr, lr_idx,~] = unique([data.lr]','stable');

steps = [data.step]';

figure(plot_num)
set(gcf,'position',([100,100,1400,600]),'color','w')
hold on
box on
grid on
plot(steps, [data.tr_ssim]','.-b');
plot(steps, [data.te_ssim]','.-g');
stem(steps(lr_idx(2:end)), ones(1,numel(lr_idx(2:end))),'.r');
set(gca,'fontweight','bold')
xlim([steps(1) steps(end)]);
ylim([0 1]);
xlabel('Training Steps', 'fontweight','bold');
ylabel('SSIM Value', 'fontweight','bold');
title('Training and Testing SSIM Results', 'fontweight','bold');
legend('Training', 'Testing', 'Learning Rate Transition', 'location', 'southoutside', 'orientation','horizontal');
plot_num = plot_num + 1;
ax = gca;
ax.Position = [0.05 0.15 0.93 0.80];

max_nrmse = max(max([data.tr_nrmse]',[data.te_nrmse]'));
max_nrmse = (ceil(max_nrmse*10)/10);
figure(plot_num)
set(gcf,'position',([100,100,1400,600]),'color','w')
hold on
box on
grid on
plot(steps, [data.tr_nrmse]','.-b');
plot(steps, [data.te_nrmse]','.-g');
stem(steps(lr_idx(2:end)), max_nrmse*ones(1,numel(lr_idx(2:end))),'.r');
set(gca, 'fontweight','bold')
xlim([steps(1) steps(end)]);
%xticks([1:length(E_name{1})]);
%xticklabels(E_name{1});
ylim([0 max_nrmse]);
xlabel('Training Steps', 'fontweight','bold');
ylabel('NRMSE Value', 'fontweight','bold');
title('Training and Testing NRMSE Results', 'fontweight','bold');
legend('Training', 'Testing', 'Learning Rate Transition', 'location', 'southoutside', 'orientation','horizontal');

ax = gca;
ax.Position = [0.05 0.15 0.93 0.80];
plot_num = plot_num + 1;


max_nmae = max(max([data.tr_nmae]',[data.te_nmae]'));
max_nmae = (ceil(max_nmae*10)/10);
figure(plot_num)
set(gcf,'position',([100,100,1400,600]),'color','w')
hold on
box on
grid on
plot(steps, [data.tr_nmae]','.-b');
plot(steps, [data.te_nmae]','.-g');
stem(steps(lr_idx(2:end)), max_nmae*ones(1,numel(lr_idx(2:end))),'.r');
set(gca, 'fontweight','bold')
xlim([steps(1) steps(end)]);
%xticks([1:length(E_name{1})]);
%xticklabels(E_name{1});
ylim([0 max_nmae]);
xlabel('Training Steps', 'fontweight','bold');
ylabel('NMAE Value', 'fontweight','bold');
title('Training and Testing NMAE Results', 'fontweight','bold');
legend('Training', 'Testing', 'Learning Rate Transition', 'location', 'southoutside', 'orientation','horizontal');

ax = gca;
ax.Position = [0.05 0.15 0.93 0.80];
plot_num = plot_num + 1;

%% get the first image in the training saves
file_filter = {'*.png','PNG Image File';'*.*','All Files' };

% read in the log file recorded during the learning process
startpath = 'D:\IUPUI\PhD\Results';
[img_file, img_path] = uigetfile(file_filter, 'Select First Training Image File', log_path);
if(img_path == 0)
    return;
end

%% cyce through the images calculated during training

% get the listing of the images
listing = dir(img_path);
listing = listing(3:end);
list_size = size(listing,1);
fprintf('Number of images: %d\n',list_size);

img_root_name = strsplit(img_file,'000500.png');
img_root_name = img_root_name{1};
im_rt_length = length(img_root_name);
img_step = [];
img_name = {};

index = 1;
for idx=1:list_size
    
    if(strncmp(listing(idx).name, img_root_name, im_rt_length))
        tmp = strsplit(listing(idx).name,{img_root_name,'.png'});
        img_step(index,1) = str2double(tmp{2});
        
        img_name{index,1} = listing(idx).name;
        index = index + 1;
    
    end
    
end

[img_step, img_idx] = sort(img_step);
img_name = img_name(img_idx,1);

img_path = listing(1).folder;

fig = figure(plot_num);
set(gcf,'position',([100,100,400,350]),'color','w')
for idx=1:length(img_step)
    
    tmp_img = imread(fullfile(img_path,img_name{idx}));
    
    image(tmp_img);
    colormap(gray(256));
    axis off
    title(strcat('DfD-Net Training Step',32,num2str(img_step(idx),'%06d')),'fontweight','bold');
    ax = gca;
    ax.Position = [0.00 0.0 1 0.94];
    
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);    
end

bp = 1;

gif_file_name = fullfile(img_path,strcat(img_root_name,'animated.gif'));
delay = 0.25*ones(length(img_step),1);
loop_count = 0;
create_animated_gif(gif_file_name, delay, im, loop_count)



