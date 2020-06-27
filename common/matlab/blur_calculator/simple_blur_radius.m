format long g
format compact
%clc
close all


% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% Setup the lens and camera parameters
px_size = 0.0048;       % pixel size (mm)
fl = 9.6;            % mm
f_num = 3.7;        % unitless
d_o = 1.0*1000*1000;     % mm

limits =[100, 1000];  % m



%% run the parameters through the blurCalc function and plot the results

c_lim = 1*px_size;
Dn = ((d_o)*(fl*fl)/(fl*fl+c_lim*f_num*(d_o-fl)))/1000;
Df = ((d_o)*(fl*fl)/(fl*fl-c_lim*f_num*(d_o-fl)))/1000;

DOF = Df-Dn;

[S_range, CoC, CoC_max] = blurCalc(f_num, fl, d_o, limits*1000);
px = ceil(CoC/px_size);
y_axis_ticks = [0:px_size:(CoC_max+px_size)];
y_axis_labels = num2str(y_axis_ticks'/px_size);

%% Plot the main blur radius curve
save_location = 'D:\IUPUI\PhD\Images\Optics';
lw = 1.0;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
plot(S_range/1000, CoC, 'LineWidth', lw, 'LineStyle', '-', 'Marker', '.', 'Color', 'b', 'MarkerSize', 9);
plot(S_range/1000, px*px_size, 'LineWidth', lw, 'LineStyle', '-', 'Marker', 'none', 'Color', 'k', 'MarkerSize', 9);
plot([0, S_range(end)/1000], [CoC_max, CoC_max], 'LineWidth', lw, 'LineStyle', '--', 'Marker', 'none', 'Color', 'g', 'MarkerSize', 9);
stem([Dn, Df],[100, 100], 'LineWidth', lw, 'LineStyle', '-', 'Marker', 'none', 'Color', 'r', 'MarkerSize', 9);
set(gca,'FontSize', 13, 'fontweight','bold');

% X-Axis
xlim(limits);
xlabel('Distance From Lens (m)', 'fontweight','bold','FontSize', 13);

% Y-Axis
ylim([0, ceil(CoC_max/px_size)*px_size]);
yticks(y_axis_ticks);
yticklabels(y_axis_labels);
ylabel('Blur Radius (pixels)', 'fontweight','bold','FontSize', 13);

title('Object Distance vs. Radius of Blur', 'fontweight','bold','FontSize', 16);
legend('Blur Radius','Quantized Blur Radius',strcat('CoC_{max} =', 32, num2str(CoC_max),'mm'),'Depth of Field', 'location', 'southeast');
ax = gca;
ax.Position = [0.05 0.11 0.93 0.83];

%print(plot_num, '-dpng', fullfile(save_location,'blur_radius2.png'));

plot_num = plot_num + 1;


%% quantize the results for a given pixel size

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
%plot(S_range/1000, CoC/px_size,'.-b');
plot(S_range/1000, px,'-k');
plot([0, S_range(end)/1000], ceil([CoC_max/px_size, CoC_max/px_size]),'-r');

set(gca, 'xlim', limits, 'fontweight','bold');
set(gca, 'ylim', [0, px(end) + 2], 'fontweight','bold');

xlabel('Distance From Lens (m)', 'fontweight','bold');
ylabel('Blur Radius (pixels)', 'fontweight','bold');
title('Object Distance vs. Radius of Blur', 'fontweight','bold');

legend(strcat('Quantized Blur Radius px_{size}:',32, num2str(px_size*1000), '{\mu}m'),strcat('R_{max} =', 32, num2str(ceil(CoC_max/px_size))), 'location', 'southeast');

ax = gca;
ax.Position = [0.05 0.11 0.93 0.83];

plot_num = plot_num + 1;