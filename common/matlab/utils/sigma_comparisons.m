format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

max_classes = 256;
max_sigma = 2.5;

%% Elliot sigma like behavior


x=[max_sigma/max_classes:max_sigma/max_classes:max_sigma];
y1 = max_sigma*2*x./(1+2*abs(x));
y2 = max_sigma*0.2*x./(max_sigma-0.8*abs(x));

figure(plot_num);
hold on
plot([0:max_classes-1],y1, '.-b');
plot([0:max_classes-1],y2, '.-g');

