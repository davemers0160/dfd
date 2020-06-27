format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% load the mat file


start_path = 'D:\IUPUI\DFD\dfd_dnn_pso\pso_results\dfd_dnn';

mat_file_filter = {'*.mat','Mat Files';'*.*','All Files' };
[mat_save_file, mat_save_path] = uigetfile(mat_file_filter, 'Select Mat File', start_path);
if(mat_save_path == 0)
    return;
end

load(fullfile(mat_save_path,mat_save_file));

commandwindow;


%% create the patches

con_num = numel(X(1,1).con(:,1));
num_itr = itr;
x_off = 4.5;
y_off = 4.5;

c = cell(num_itr, con_num);

for idx=1:num_itr-1
    for jdx=1:con_num
    
        c{idx,jdx} = make_conv_patch(G(idx).con(jdx,:), x_off+10*(idx-1), y_off);
    end
end

%% plot the two convs for G

% figure(plot_num)
% set(gcf,'position',([50,50,1200,700]),'color','w')
% 
% subplot(2,1,1);
% hold on
% box on
% grid on
% 
% for idx=1:num_itr
% 
%     for jdx=1:con_num
%         
%         patch('Faces',[1 2 3 4], 'Vertices',c{idx,jdx},'FaceColor','blue','FaceAlpha',0.3, 'LineStyle','-','EdgeColor','b', 'LineWidth', 1);
% 
%     end
% end
% 
% 
% xlim([0,10*num_itr + 5]);
% xticks([4.5:10:10*num_itr]);
% xticklabels(num2str([1:1:num_itr]'));
% 
% ylim([0, 10*con_num]);
% yticks([4.5:10:(10*(con_num-1)+4.5)]);
% %yticklabels({'con 2', 'con 1'});
% 
% 
% subplot(2,1,2);
% 
% plot_num = plot_num + 1;

%%

N = size(P,1);
num_itr = itr;
X_con = [];
X_cond = [];
X_cont = [];
X_act = [];
X_bn = [];
X_crop = [];

con_index = find(strcmp('con',net_description.net_structure));
% cond_index = find(strcmp('cond',net_description.net_structure));
% cont_index = find(strcmp('cont2u',net_description.net_structure));

con_num = numel(X(1,1).con(:,1));
% cond_num = numel(cond_index);
% cont_num = numel(cont_index);
act_num = numel(X(1,1).act(:,1));
bn_num = numel(X(1,1).bn(:,1));

for idx=1:N
   
    for jdx=1:num_itr

        for kdx = 1:con_num
            X_con(idx,jdx,kdx) = X(idx,jdx).con(kdx,1);
        end
        
        for kdx=1:act_num
            X_act(idx,jdx,kdx) = X(idx,jdx).act(kdx);
        end
        
        for kdx=1:bn_num
            X_bn(idx,jdx,kdx) = X(idx,jdx).bn(kdx);
        end
        
        X_crop(idx,jdx) = X(idx,jdx).crop_size;
        
%         for kdx = 1:cont_num
%             X_con(idx,jdx,kdx) = X(idx,jdx).con(kdx,1);
%         end
%         
%         for kdx = 1:cont_num
%             X_con(idx,jdx,kdx) = X(idx,jdx).con(kdx,1);
%         end
        
    end
    
end

%% plots


for idx=2:con_num
    
    figure(plot_num);
    set(gcf,'position',([50,50,1200,600]),'color','w')
    box on;
    grid on;
    surf(X_con(:,:,idx));

    set(gca,'fontweight','bold', 'fontsize',13);

    xlim([0,num_itr+1]);
    xticks([0:2:num_itr]);
    %xticklabels(num2str([1:1:num_itr]'));
    xtickangle(35);
    xlabel('Iterations', 'fontweight','bold', 'fontsize',13);

    ylim([0, N+1]);
    yticks([0:1:N]);
    %yticklabels({'con 2', 'con 1'});
    ylabel('Particle Number', 'fontweight','bold', 'fontsize',13);


    zlim([0 525]);
    zticks([0:25:525]);
    ztickformat('%2.1f');
    zlabel('Number of Filters', 'fontweight','bold', 'fontsize',13);
    
    title(strcat(net_description.net_structure(con_index(idx)), num2str(con_index(idx),': %02d')), 'fontweight','bold', 'fontsize',13);

    view(45,30);
    ax = gca;
    ax.Position = [0.08 0.1 0.9 0.86];
    plot_num = plot_num + 1;

end

%%

for idx=1:act_num
    
    figure(plot_num);
    set(gcf,'position',([50,50,1200,600]),'color','w')
    box on;
    grid on;
    surf(X_act(:,:,idx));

    set(gca,'fontweight','bold', 'fontsize',13);

    xlim([0,num_itr+1]);
    xticks([0:2:num_itr]);
    %xticklabels(num2str([1:1:num_itr]'));
    xtickangle(35);
    xlabel('Iterations', 'fontweight','bold', 'fontsize',13);

    ylim([0, N+1]);
    yticks([0:1:N]);
    %yticklabels({'con 2', 'con 1'});
    ylabel('Particle Number', 'fontweight','bold', 'fontsize',13);


    zlim([0 7]);
    zticks([0:1:7]);
    ztickformat('%2.1f');
    zlabel('Number of Filters', 'fontweight','bold', 'fontsize',13);

    view(45,30);
    ax = gca;
    ax.Position = [0.08 0.1 0.9 0.88];
    plot_num = plot_num + 1;
    
end


%%

for idx=1:bn_num
    
    figure(plot_num);
    set(gcf,'position',([50,50,1200,600]),'color','w')
    box on;
    grid on;
    surf(X_bn(:,:,idx));

    set(gca,'fontweight','bold', 'fontsize',13);

    xlim([0,num_itr+1]);
    xticks([0:2:num_itr]);
    %xticklabels(num2str([1:1:num_itr]'));
    xtickangle(35);
    xlabel('Iterations', 'fontweight','bold', 'fontsize',13);

    ylim([0, N+1]);
    yticks([0:1:N]);
    %yticklabels({'con 2', 'con 1'});
    ylabel('Particle Number', 'fontweight','bold', 'fontsize',13);


    zlim([0 2]);
    zticks([0:1:2]);
    ztickformat('%2.1f');
    zlabel('Number of Filters', 'fontweight','bold', 'fontsize',13);

    view(45,30);
    ax = gca;
    ax.Position = [0.08 0.1 0.9 0.88];
    plot_num = plot_num + 1;
    
end

%%

figure(plot_num);
set(gcf,'position',([50,50,1200,600]),'color','w')
box on;
grid on;
surf(X_crop(:,:));

set(gca,'fontweight','bold', 'fontsize',13);

xlim([0,num_itr+1]);
xticks([0:2:num_itr]);
%xticklabels(num2str([1:1:num_itr]'));
xtickangle(35);
xlabel('Iterations', 'fontweight','bold', 'fontsize',13);

ylim([0, N+1]);
yticks([0:1:N]);
%yticklabels({'con 2', 'con 1'});
ylabel('Particle Number', 'fontweight','bold', 'fontsize',13);


zlim([2 16]);
zticks([2:2:16]);
ztickformat('%2.1f');
zlabel('Crop Size (Pixels)', 'fontweight','bold', 'fontsize',13);

view(45,30);
ax = gca;
ax.Position = [0.08 0.1 0.9 0.88];
plot_num = plot_num + 1;
    
%% Plot the path of a single parameter within a particle

particle_num = 2;
param_index = 9;
particle_name = strcat('Conv Layer',32,num2str(con_index(param_index),': %02d'));

p = X_con(particle_num,:,param_index);
pb = [P(particle_num,:).con];
pb = pb(param_index,1:3:end);

gb = [G.con];
gb = gb(param_index,1:3:end);

y_plt_step = 10;

x_min = 25;
x_max = num_itr;
y_min = floor(min(p(x_min:end))/y_plt_step)*y_plt_step;
y_max = ceil(max(p(x_min:end))/y_plt_step)*y_plt_step;

figure(plot_num)
set(gcf,'position',([50,50,1200,600]),'color','w')
hold on
grid on
box on

plot([1:num_itr], p, 'LineStyle','-', 'LineWidth', 1, 'Color','b', 'Marker', '.');

set(gca,'fontweight','bold', 'fontsize',13);

xlim([x_min x_max]);
xticks([x_min:2:x_max]);
xtickangle(90);
xlabel('Iteration', 'fontweight', 'bold', 'FontSize', 13);
set(gca,'XMinorTick','on', 'XMinorGrid','on');

ylim([y_min y_max]);
yticks([y_min:y_plt_step:y_max]);
ylabel('Value', 'fontweight', 'bold', 'FontSize', 13);
set(gca,'YMinorTick','on', 'YMinorGrid','on');

title(strcat('Particle',32,num2str(particle_num),32,'-',32,particle_name,32,'Convergence Track'),'fontweight','bold','FontSize', 16)

ax = gca;
ax.Position = [0.06 0.12 0.92 0.8];
ax.XAxis.MinorTickValues = x_min:1:x_max;
ax.YAxis.MinorTickValues = y_min:1:y_max;

plot_num = plot_num + 1;

%%
figure(plot_num);
set(gcf,'position',([50,50,1200,600]),'color','w')
box on;
grid on;

plot(g_best,'.-b')

plot_num = plot_num + 1;




