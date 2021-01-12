format long g
format compact
%clc
close all


% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% Setup the lens and camera parameters
f = 9.6;            % mm
f_num = 3.7;        % unitless
focus_dist = 0.5;   % m
px_size = 0.0048;        % mm
limits =[0.1, 10];  % m

step = 10;
df = focus_dist*1000;
CoC_max = (f*f/(f_num*(df-f)));

l = length(df:step:limits(2)*1000);

% D = -0.01*ones(l,l);
idx=1;
jdx=1;

% for D1=df:step:limits(2)*1000
%     jdx=1;
%     for D2=D1:step:limits(2)*1000
%         
%         D(idx,jdx) = ((df*f*f)/(f_num*(df-f)))*((1/D1) - (1/D2));
%         jdx = jdx + 1;
%     end
%     idx = idx + 1;
% end
D=[];

D1=df:step:limits(2)*1000;

D = ((df*f*f)/(f_num*(df-f)))*(1./D1);

range_step = 1;

[S_range, CoC, CoC_max] = blurCalc(f_num, f, focus_dist, limits, range_step);

%%
figure(plot_num)
set(gcf,'position',([100,100,800,600]),'color','w')
hold on
box on
grid on

plot(D1/1000,D/px_size,'.-b')
plot(S_range/1000, CoC/px_size,'.-g');

set(gca, 'xlim', [0,6],'fontweight','bold')
set(gca, 'ylim', [1,ceil(max(D/px_size))],'fontweight','bold')

% surf(D)
% shading interp
% colormap(jet(1000))
% 
% z_axis_ticks = [px:px:(CoC_max+px)];
% z_axis_labels = num2str(z_axis_ticks'/px);
% set(gca, 'zlim',[px,max(D(:))], 'zTickLabel', z_axis_labels, 'ZTick',z_axis_ticks);
% % set(gca,'Xlim', limits);
% % set(gca,'Ylim', limits);

plot_num = plot_num + 1;