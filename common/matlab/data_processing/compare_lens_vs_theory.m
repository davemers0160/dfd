format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;
global pixel px_size range

%% Select the lens info file that describes the number of blurred pixels

file_filter = {'*.txt','Text Files';'*.*','All Files' };

startpath = 'D:\IUPUI\DfD\dfd_dnn_rw';
[data_file, data_path] = uigetfile(file_filter, 'Select Configuration File', startpath);
if(data_path == 0)
    return;
end

commandwindow;

%% process the lens info

div = 7;
data_params = parse_input_parameters(fullfile(data_path, data_file));

% row 1 is the pixel blur size
range = str2double({data_params{1}{2:end}})*1000;
limits =[0.1, max(range)];  % m
px_size = 0.0048;                       % pixel size (mm)
 
data_params(1) = [];

num_steps = numel(data_params);

lens_step = cell(num_steps,1);
lens_pixel = cell(num_steps,1);

for idx=1:num_steps
    
    lens_step{idx,1} = data_params{idx}{1};
    lens_pixel{idx,1} = str2double({data_params{idx}{2:end}});
    
end


%% match the rw lens value to the theoretical

%f_num, focal length, d_o
x_lim = [1.0,5.0; 9.88,9.88; 14000,16000];
v_max = [-0.02, 0.02; -0.0, 0.0; -5, 5];
itr_max = 600;
N = 5000;
c1 = 2.1;
c2 = 2.1;

c_lim = 1*px_size;
%d_o = 1.2*1000;                         % mm

f_num = zeros(num_steps,1);
fl = zeros(num_steps,1);
d_o = zeros(num_steps,1);

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')

l_start = 1;
l_stop = 2;     % num_steps

for idx=l_start:l_stop

    fprintf('\n');
    fprintf('Voltage Step, Focal Distance ${d_{o}}$, ${F_{num}}$, Focal Length, Error\n');
    
    pixel = lens_pixel{idx,1};   
    tmp_do = range(pixel<=2);  
    
%     if(isempty(tmp_do))
%         x_lim(3,:) = [40, 500000];
%         
%     elseif(numel(tmp_do) == 1)
%         if(tmp_do == range(1))
%             x_lim(3,:) = [40, tmp_do];
%         elseif(tmp_do == range(end))
%             x_lim(3,:) = [tmp_do, 500000];
%         end
%     else
%         if(tmp_do(1) == tmp_do(2))
%             x_lim(3,:) = [40, 500000];
%         else
%             x_lim(3,:) = [tmp_do(1), tmp_do(2)];
%         end
%     end
    
    for jdx=1:4

        [x, v, g, p, pso_stats, itr_cnt] = PSO(@get_coc, N, itr_max, v_max, x_lim, 1.7, c1, c2, 'constrict');

        [err] = get_coc(g(:,itr_cnt)');

        f_num(idx,1) = g(1, itr_cnt);
        fl(idx,1) = g(2,itr_cnt);
        d_o(idx,1) = g(3,itr_cnt);

        fprintf('%02d, %s, %2.4f, %2.4f, %2.4f, %2.4f\n', jdx, lens_step{idx,1}, d_o(idx,1), f_num(idx,1), fl(idx,1), err);
        
        
        Dn = ((d_o(idx,1))*(fl(idx,1)*fl(idx,1))/(fl(idx,1)*fl(idx,1)+c_lim*f_num(idx,1)*(d_o(idx,1)-fl(idx,1))))/1000;
        Df = ((d_o(idx,1))*(fl(idx,1)*fl(idx,1))/(fl(idx,1)*fl(idx,1)-c_lim*f_num(idx,1)*(d_o(idx,1)-fl(idx,1))))/1000;

        DOF = Df-Dn;

        [S_range, CoC, CoC_max] = blurCalc(f_num(idx,1), fl(idx,1), d_o(idx,1), limits);
        px = ceil(CoC/px_size);
        
        plot(S_range/1000, CoC, 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','b', 'MarkerSize', 8);

        hold on;
        box on;
        grid on;
        plot(S_range/1000, px*px_size, 'Color','k', 'LineStyle','-', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','k', 'MarkerSize', 10);
        % plot([0, S_range(end)/1000], [CoC_max, CoC_max],'-g');
        % stem([Dn, Df],[CoC_max, CoC_max],'.r');

        %plot(lens_pixel{step_index,1}/1000, pixel(1:numel(lens_pixel{step_index,1}))*px_size, 'Color','g', 'LineStyle','-', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','g', 'MarkerSize', 10);
        plot(range/1000, pixel*px_size, 'Color','g', 'LineStyle','--', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','g', 'MarkerSize', 10);

        set(gca,'FontSize', 13, 'fontweight','bold');

        % X-Axis
        xlim(limits/1000);
        xticks(limits(1)/1000:0.2:limits(2)/1000);
        xtickformat('%2.1f');
        xlabel('Distance From Lens (m)', 'fontweight','bold','FontSize', 13);

        % Y-Axis
        %y_axis_ticks = [0:px_size:(CoC_max+px_size)];
        y_axis_ticks = [0:px_size:(max(pixel)*px_size+px_size)];
        y_axis_labels = num2str(y_axis_ticks'/px_size);

        %ylim([0, ceil(CoC_max/px_size)*px_size]);
        ylim([0, (max(pixel)+1)*px_size]);
        yticks(y_axis_ticks);
        yticklabels(y_axis_labels);
        ylabel('Blur Radius (pixels)', 'fontweight','bold','FontSize', 13);

        title(strcat('Object Distance vs. Radius of Blur - Voltage Step:',32,lens_step{idx,1}), 'fontweight','bold', 'FontSize', 16);
        legend('Theoretical Blur Radius','Quantized Blur Radius','Measured Lens Blur Radius', 'location', 'southoutside', 'orientation','horizontal');

        str = {sprintf('f-number = %2.2f', f_num(idx,1)),...
               sprintf('focus distance = %2.3fm', d_o(idx,1)/1000)};
        delete(findall(gcf,'Tag','info'));
        annotation('textbox',[0.76,0.81,0.2,0.09],'String',str,'FitBoxToText','off','fontweight','bold', 'FontSize', 12, 'BackGroundColor','w','Tag','info');
        ax = gca;
        ax.Position = [0.05 0.16 0.93 0.77];
        drawnow;
        hold off;
    end
end

plot_num = plot_num + 1;

return;

%% Setup the lens and camera parameters

% Dn = ((d_o)*(fl*fl)/(fl*fl+c_lim*f_num*(d_o-fl)))/1000;
% Df = ((d_o)*(fl*fl)/(fl*fl-c_lim*f_num*(d_o-fl)))/1000;
% 
% DOF = Df-Dn;

fl2 = 9.88;

% 141
% f_num2 = 3.9550; % 3.9550, 3.7360;
% d_o2 = 499.8000; % 499.8000, 526.2000 

% 129
f_num2 = 2.8940;
d_o2 = 126658.0000; 

[S_range, CoC, CoC_max] = blurCalc(f_num2, fl2, d_o2, limits);
px = ceil(CoC/px_size);

l_start = 1;
pixel = lens_pixel{l_start,1}; 

%% Plot the main blur radius curve
save_location = 'D:\IUPUI\PhD\Images\Camera';

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
plot(S_range/1000, CoC, 'Color','b', 'LineStyle','-', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','b', 'MarkerSize', 8);
plot(S_range/1000, px*px_size, 'Color','k', 'LineStyle','-', 'LineWidth',1, 'Marker','none', 'MarkerFaceColor','k', 'MarkerSize', 10);
% plot([0, S_range(end)/1000], [CoC_max, CoC_max],'-g');
% stem([Dn, Df],[CoC_max, CoC_max],'.r');

%plot(lens_pixel{step_index,1}/1000, pixel(1:numel(lens_pixel{step_index,1}))*px_size, 'Color','g', 'LineStyle','-', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','g', 'MarkerSize', 10);
plot(range/1000, pixel*px_size, 'Color','g', 'LineStyle','--', 'LineWidth',1, 'Marker','.', 'MarkerFaceColor','g', 'MarkerSize', 10);

set(gca,'FontSize', 13, 'fontweight','bold');

% X-Axis
xlim(limits/1000);
xticks(limits(1)/1000:0.2:limits(2)/1000);
xtickformat('%2.1f');
xlabel('Distance From Lens (m)', 'fontweight','bold','FontSize', 13);

% Y-Axis
%y_axis_ticks = [0:px_size:(CoC_max+px_size)];
y_axis_ticks = [0:px_size:(max(pixel)*px_size+px_size)];
y_axis_labels = num2str(y_axis_ticks'/px_size);

%ylim([0, ceil(CoC_max/px_size)*px_size]);
ylim([0, (max(pixel)+1)*px_size]);
yticks(y_axis_ticks);
yticklabels(y_axis_labels);
ylabel('Blur Radius (pixels)', 'fontweight','bold','FontSize', 13);

title(strcat('Object Distance vs. Radius of Blur - Voltage Step:',32,lens_step{l_start,1}), 'fontweight','bold', 'FontSize', 16);
legend('Theoretical Blur Radius','Quantized Blur Radius','Lens Blur Radius', 'location', 'southoutside', 'orientation','horizontal');

str = {sprintf('f-number = %2.2f', f_num2),...
       sprintf('focal distance = %2.3fm', d_o2/1000)};
delete(findall(gcf,'Tag','info'));
annotation('textbox',[0.76,0.81,0.2,0.09],'String',str,'FitBoxToText','off','fontweight','bold', 'FontSize', 12, 'BackGroundColor','w','Tag','info');
ax = gca;
ax.Position = [0.05 0.16 0.93 0.77];

plot_num = plot_num + 1;

return

%% test things

% x1 = [x_lim(1,1):0.005:x_lim(1,2)];
% x3 = [x_lim(3,1):0.2:x_lim(3,2)];
pixel = lens_pixel{2,1};

%x1 = [2.3:0.1:3.4];    1472725.7
x1 = [0.20:0.010:5.0];
x3 = [100.0:0.005:600];

x2 = ones(1,numel(x3))*9.88;

err = [];
tic;
for idx=1:numel(x1)
    x = cat(1, x1(idx)*ones(1,numel(x3)),x2,x3);
    [err(idx,:)] = get_coc(x(:,:))';
end
toc


me = 51;
max_err = max(err(:));
err(err>me) = me;

figure(plot_num)
surf(x3, x1, err)
shading interp
xlabel('d_o');
ylabel('f_{num}');
zlim([0,me-1]);
% zlim([0,20]);
colormap(jet(me));
view(90,90);

fprintf('f-num: %2.3f - %2.3f\n', x1(1), x1(end));
fprintf('Range: %06d - %06d\n', x3(1), x3(end));
fprintf('Min: %d\n', min(err(:)));
fprintf('Max: %d\n\n', max_err);
plot_num = plot_num + 1;

%% ----------------------------

% function [err] = get_coc(x)
%     global pixel px_size range
% 
%     
%     f_num = x(1);
%     f = x(2);
%     d_o = x(3);
%     
%     %     px_size = 0.0048;    % pixel size (mm)
% % %     pixel = [7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9];
% %     
%     
%     % 135
% %     d_o = 0.93*1000;     % mm
% %     d_near = [220, 248, 308, 346, 414, 482, 754];
% %     d_far = [1232, 1910, 3589];
% %     em = [2 1 1 1 1 1 1 4 10 12];
% %     pixel = pixel(1:10);
%     
%     % 129
% %     d_o = 5.5*1000;     % mm
% %     d_near = [506,610,740,926,1168,1744,5110];
% %     d_far = [];  
%     
%     dn = (range <= d_o);
%     
%     d_near = range(dn);
%     d_far = range(~dn);
%     
% %     em = [2 2 2 1 1 1 1];
% 
% 
%     % 143
% %     d_o = 0.262*1000;     % mm
% %     d_near = [110, 134, 152, 168, 190, 206, 250];
% %     d_far = [296, 398, 470, 548, 672, 914, 1260, 1934, 2970];  
% %     em1 = [2 2 2 1 1 1 1];
% %     em2 = [1 1 1 1 1 1 1 1 1];
% %     em = cat(2,em1,em2);
% %     pixel = pixel(1:16);
%     
% 
%             
%     tl = (d_o*f*f)/(f_num*(d_o-f));
%     
%     coc_far = tl*((1/d_o)-(1./d_far));
%     coc_near = tl*((1./d_near) - (1/d_o));
%     
%     CoC = [coc_near coc_far];  
%     
%     
%     px = (CoC/px_size);
%     px = ceil(px);
%     
% %     err = (pixel - px).*(pixel - px);
% %     err = mean(err);
% %     err = sum(abs(pixel - px).*em);
%     err = sum(abs(pixel - px));
%     
% end


