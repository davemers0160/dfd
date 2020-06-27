format long g
format compact
clc
clearvars
close all


% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% Start looking at the difference between the two methods for a giving dataset
save_dir = 'D:\IUPUI\PhD\Images\gc_dnn_comp\';
gc_results_dir = 'D:\IUPUI\PhD\Results\dfd_gc';
dnn_results_dir = 'D:\IUPUI\PhD\Results\dfd_dnn';

gc_results_dir = uigetdir(gc_results_dir,'Select Graph Cuts Result Folder');
if(gc_results_dir == 0)
    return;
end

% load the graph cuts results
load(fullfile(gc_results_dir,'dfd_gc_results.mat'), 'data', 'data_name');

%rename the variables
gc_data_i = data;
gc_data_name_i = data_name;
% gc_data_sort = data_sort;
% gc_data_name_sort = name_sort;
% gc_s_i = s_i;


dnn_results_dir = uigetdir(dnn_results_dir,'Select DFD-Net Result Folder');
if(dnn_results_dir == 0)
    return;
end
% load the dnn results
load(fullfile(dnn_results_dir,'dfd_dnn_test_results.mat'), 'data', 'data_name');

%rename the variables
dnn_data = data;
dnn_data_name = data_name;
% dnn_data_sort = data_sort;
% dnn_data_name_sort = name_sort;
% dnn_s_i = s_i;

clear data data_name s_i data_sort name_sort

%% Since the dnn results are smaller than the gc results we can only compare the dnn test dataset

%hard coded right now
gc_idx = [19:36, 91:108, 369:386];

gc_data = gc_data_i(gc_idx,:);
gc_data_name = gc_data_name_i(gc_idx,:);

%% plot the results for each illumination value

iNMAE = 4;
iNRMSE = 5;
iSSIM = 6;
irun = 7;

% for illum=1:3
%     gc_plot_I = gc_data(gc_data(:,1) == illum,:);
%     gc_I_name = gc_data_name(gc_data(:,1) == illum,:);
%     
%     dnn_plot_I = dnn_data(dnn_data(:,1) == illum,:);
%     dnn_I_name = dnn_data_name(dnn_data(:,1) == illum,:);
%     
%     gc_plot_E = {};
%     gc_E_name = {};
%     
%     for idx=1:3
%         gc_plot_E{idx} = gc_plot_I(gc_plot_I(:,2) == idx,:);
%         dnn_plot_E{idx} = dnn_plot_I(dnn_plot_I(:,2) == idx,:);
%         
%         gc_E_name{idx} = gc_I_name(gc_plot_I(:,2) == idx,:);
%     end
%     
%     fig = figure(plot_num);
%     set(gcf,'position',([50,50,1400,700]),'color','w')
%     set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
%     
%     nrmse_max = max(max(gc_data(:,iNRMSE)),max(dnn_data(:,iNRMSE)));
%     subplot(3,1,1)
%     hold on
%     box on
%     grid on
% 
%     set(gca,'fontweight','bold')
%     yyaxis left
%     plot(gc_plot_E{1}(:,iNRMSE))
%     plot(gc_plot_E{2}(:,iNRMSE))
%     plot(gc_plot_E{3}(:,iNRMSE))
%     ylim([0 ceil(nrmse_max)])
%     ylabel('Graph Cuts','fontweight','bold')
%     title(strcat('NRMSE Comparision - Illumination',32,num2str(illum)),'fontweight','bold')
%     
%     yyaxis right
%     plot(dnn_plot_E{1}(:,iNRMSE))
%     plot(dnn_plot_E{2}(:,iNRMSE))
%     plot(dnn_plot_E{3}(:,iNRMSE))
%     ylim([0.0 ceil(nrmse_max)])
%     ylabel('DFD-Net','fontweight','bold')
% 
%     xlim([1 length(gc_E_name{1})]);
%     xticks([1:length(gc_E_name{1})]);
%     %xticklabels(gc_E_name{1});
%     xticklabels([]);
%     xtickangle(90);
%     ax = gca;
%     ax.Position = [0.05 0.72 0.9 0.24];
%     
%     nmae_max = max(max(gc_data(:,iNMAE)),max(dnn_data(:,iNMAE)));
%     subplot(3,1,2)
%     hold on
%     box on
%     grid on
% 
%     set(gca,'fontweight','bold')
%     yyaxis left
%     plot(gc_plot_E{1}(:,iNMAE))
%     plot(gc_plot_E{2}(:,iNMAE))
%     plot(gc_plot_E{3}(:,iNMAE))
%     ylim([0 (nmae_max+.1)])
%     ylabel('Graph Cuts','fontweight','bold')
%     title(strcat('NMAE Comparision - Illumination',32,num2str(illum)),'fontweight','bold')
% 
%     yyaxis right
%     plot(dnn_plot_E{1}(:,iNMAE))
%     plot(dnn_plot_E{2}(:,iNMAE))
%     plot(dnn_plot_E{3}(:,iNMAE))
%     ylim([0.0 (nmae_max+.1)])
%     ylabel('DFD-Net','fontweight','bold')
% 
%     xlim([1 length(gc_E_name{1})]);
%     xticks([1:length(gc_E_name{1})]);
%     %xticklabels(gc_E_name{1});
%     xticklabels([]);
%     xtickangle(90);
%     ax = gca;
%     ax.Position = [0.05 0.435 0.90 0.24];
% 
%     subplot(3,1,3)
%     hold on
%     box on
%     grid on
% 
%     set(gca,'fontweight','bold')
%     yyaxis left
%     plot(gc_plot_E{1}(:,iSSIM))
%     plot(gc_plot_E{2}(:,iSSIM))
%     plot(gc_plot_E{3}(:,iSSIM))
%     ylim([0 1])
%     ylabel('Graph Cuts','fontweight','bold')
%     title(strcat('SSIM Comparision - Illumination',32,num2str(illum)),'fontweight','bold')
% 
%     yyaxis right
%     plot(dnn_plot_E{1}(:,iSSIM))
%     plot(dnn_plot_E{2}(:,iSSIM))
%     plot(dnn_plot_E{3}(:,iSSIM))
%     ylim([0.0 1.0])
%     ylabel('DFD-Net','fontweight','bold')
% 
%     xlim([1 length(gc_E_name{1})]);
%     xticks([1:length(gc_E_name{1})]);
%     xticklabels(gc_E_name{1});
%     xtickangle(90);
%     legend('Exposure 0', 'Exposure 1', 'Exposure 2', 'Exposure 0', 'Exposure 1', 'Exposure 2', 'location', 'southoutside', 'orientation', 'horizontal')
%     ax = gca;
%     ax.Position = [0.05 0.15 0.9 0.24];
% 
%     print(plot_num, '-dpng', fullfile(save_dir,strcat('DNN_GC_Comp_I',num2str(illum),'.png')));
% 
%     plot_num = plot_num + 1;       
%     
%     
% end

%% try a plot of all combined

% need to build up a name listing
% illum_name = {'I1','I2', 'I3'};
% exp_name = {'E0','E1','E2'};
% lr_name = {'L','R'};
% 
% combined_name = cell(length(dnn_data_name),1);
% for idx=1:length(dnn_data_name)
%     combined_name{idx,1} = strcat(dnn_data_name{idx}, '-', illum_name{dnn_data(idx,1)},'-',exp_name{dnn_data(idx,2)},'-',lr_name{dnn_data(idx,3)});
% end
% 
% nrmse_max = max(max(gc_data(:,iNRMSE)),max(dnn_data(:,iNRMSE)));
% fig = figure(plot_num);
% set(gcf,'position',([50,50,1400,700]),'color','w')
% set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
% hold on
% box on
% grid on
% 
% set(gca,'fontweight','bold')
% yyaxis left
% plot(gc_data(:,iNRMSE))
% plot(gc_data(:,iNMAE))
% plot(gc_data(:,iSSIM))
% ylim([0 ceil(nrmse_max)])
% ylabel('Graph Cuts','fontweight','bold')
% title('Graph Cuts vs. DFD-Net Performance Comparision','fontweight','bold')
% 
% yyaxis right
% plot(dnn_data(:,iNRMSE))
% plot(dnn_data(:,iNMAE))
% plot(dnn_data(:,iSSIM))
% ylim([0.0 ceil(nrmse_max)])
% ylabel('DFD-Net','fontweight','bold')
% 
% xlim([1 length(combined_name)]);
% xticks([1:length(combined_name)]);
% xticklabels(combined_name);
% xtickangle(90);
% legend('NRMSE', 'NMAE', 'SSIM', 'NRMSE', 'NMAE', 'SSIM', 'location', 'southoutside', 'orientation', 'horizontal')
% ax = gca;
% ax.Position = [0.05 0.22 0.9 0.73];
% 
% print(plot_num, '-dpng', fullfile(save_dir,strcat('DNN_GC_Comp_ALL.png')));
% 
% plot_num = plot_num + 1;


%% all combined, but in subplots
% 
% combined_name = cell(length(dnn_data_name),1);
% for idx=1:length(dnn_data_name)
%     combined_name{idx,1} = strcat(dnn_data_name{idx}, '-', illum_name{dnn_data(idx,1)},'-',exp_name{dnn_data(idx,2)},'-',lr_name{dnn_data(idx,3)});
% end
% nrmse_max = max(max(gc_data(:,iNRMSE)),max(dnn_data(:,iNRMSE)));
% fig = figure(plot_num);
% set(gcf,'position',([50,50,1400,700]),'color','w')
% %set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
% 
% subplot(3,1,1);
% hold on
% box on
% grid on
% set(gca,'fontweight','bold')
% plot(gc_data(:,iNRMSE),'-b');
% plot(dnn_data(:,iNRMSE),'-r');
% ylim([0 ceil(nrmse_max)]);
% ylabel('NRMSE','fontweight','bold')
% title(strcat('Graph Cut vs. DFD-Net Performance Comparision'),'fontweight','bold')
% xlim([1 length(combined_name)]);
% xticks([1:length(combined_name)]);
% xticklabels([]);
% ax = gca;
% ax.Position = [0.05 0.73 0.9 0.23];    
% 
% nmae_max = max(max(gc_data(:,iNMAE)),max(dnn_data(:,iNMAE)));
% subplot(3,1,2);
% hold on
% box on
% grid on
% set(gca,'fontweight','bold')
% plot(gc_data(:,iNMAE),'-b');
% plot(dnn_data(:,iNMAE),'-r');
% ylim([0.0 (nmae_max+.1)])
% ylabel('NMAE','fontweight','bold')
% xlim([1 length(combined_name)]);
% xticks([1:length(combined_name)]);
% xticklabels([]);
% ax = gca;
% ax.Position = [0.05 0.475 0.9 0.23];  
% 
% subplot(3,1,3);
% hold on
% box on
% grid on    
% plot(gc_data(:,iSSIM),'-b');
% plot(dnn_data(:,iSSIM),'-r')
% set(gca,'fontweight','bold')
% ylim([0.0 1.0])
% ylabel('SSIM','fontweight','bold')
% 
% xlim([1 length(combined_name)]);
% xticks([1:length(combined_name)]);
% xticklabels(combined_name);
% xtickangle(90);
% legend('Graph Cut', 'DFD-Net', 'location', 'southoutside', 'orientation', 'horizontal')
% ax = gca;
% ax.Position = [0.05 0.22 0.9 0.23];
% 
% print(plot_num, '-dpng', fullfile(save_dir,strcat('DNN_GC_Comp_ALL2.png')));
% 
% plot_num = plot_num + 1;



%% Try a similar plot but just look at one illumination value at a time


% combined_name = {};
% 
% % illum_name = {'Illum1','Illum2', 'Illum3'};
% % exp_name = {'Exp0','Exp1','Exp2'};
% % lr_name = {'Left','Right'};
% 
% illum_name = {'I1','I2', 'I3'};
% exp_name = {'E0','E1','E2'};
% lr_name = {'L','R'};
% 
% for illum=1:3
%     
%     % get the illumination specific data
%     gc_plot_I = gc_data(gc_data(:,1) == illum,:);
%     gc_I_name = gc_data_name(gc_data(:,1) == illum,:);
%     
%     dnn_plot_I = dnn_data(dnn_data(:,1) == illum,:);
%     dnn_I_name = dnn_data_name(dnn_data(:,1) == illum,:);
%     
%     combined_name = cell(length(dnn_I_name),1);
%     for idx=1:length(dnn_I_name)
%         combined_name{idx,1} = strcat(dnn_I_name{idx}, '-',exp_name{dnn_plot_I(idx,2)},'-',lr_name{dnn_plot_I(idx,3)});
%     end   
%     
%     nrmse_max = max(max(gc_plot_I(:,iNRMSE)),max(dnn_plot_I(:,iNRMSE)));
%     fig = figure(plot_num);
%     set(gcf,'position',([50,50,1400,700]),'color','w')
%     set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
%     hold on
%     box on
%     grid on
% 
%     set(gca,'fontweight','bold')
%     yyaxis left
%     plot(gc_plot_I(:,iNRMSE))
%     plot(gc_plot_I(:,iNMAE))
%     plot(gc_plot_I(:,iSSIM))
%     ylim([0 ceil(nrmse_max)])
%     ylabel('Graph Cuts','fontweight','bold')
%     title(strcat('Graph Cuts vs. DFD-Net Performance Comparision',32, '- Illumination',32,num2str(illum)),'fontweight','bold')
% 
%     yyaxis right
%     plot(dnn_plot_I(:,iNRMSE))
%     plot(dnn_plot_I(:,iNMAE))
%     plot(dnn_plot_I(:,iSSIM))
%     ylim([0.0 ceil(nrmse_max)])
%     ylabel('DFD-Net','fontweight','bold')
% 
%     xlim([1 length(combined_name)]);
%     xticks([1:length(combined_name)]);
%     xticklabels(combined_name);
%     xtickangle(90);
%     legend('NRMSE', 'NMAE', 'SSIM', 'NRMSE', 'NMAE', 'SSIM', 'location', 'southoutside', 'orientation', 'horizontal')
%     ax = gca;
%     ax.Position = [0.05 0.2 0.9 0.75];
% 
%     print(plot_num, '-dpng', fullfile(save_dir,strcat('DNN_GC_Comp2_I',num2str(illum),'.png')));
% 
%     plot_num = plot_num + 1;
% 
%     
% end

%% split them up and plot 
combined_name = {};
comp_gc_dnn = {};
illum_name = {'I1','I2', 'I3'};
exp_name = {'E0','E1','E2'};
lr_name = {'L','R'};

for illum=1:3
    
    % get the illumination specific data
    gc_plot_I = gc_data(gc_data(:,1) == illum,:);
    gc_I_name = gc_data_name(gc_data(:,1) == illum,:);
    
    dnn_plot_I = dnn_data(dnn_data(:,1) == illum,:);
    dnn_I_name = dnn_data_name(dnn_data(:,1) == illum,:);
    
    combined_name = cell(length(dnn_I_name),1);
    for idx=1:length(dnn_I_name)
        combined_name{idx,1} = strcat(dnn_I_name{idx}, '-',exp_name{dnn_plot_I(idx,2)},'-',lr_name{dnn_plot_I(idx,3)});
    end   
    
%     comp_gc_dnn{illum,1} = [gc_plot_I(:,iNRMSE), dnn_plot_I(:,iNRMSE)];
%     comp_gc_dnn{illum,2} = [gc_plot_I(:,iNMAE), dnn_plot_I(:,iNMAE)];
%     comp_gc_dnn{illum,3} = [gc_plot_I(:,iSSIM), dnn_plot_I(:,iSSIM)];
    
    nrmse_max = max(max(gc_plot_I(:,iNRMSE)),max(dnn_plot_I(:,iNRMSE)));
    fig = figure(plot_num);
    set(gcf,'position',([50,50,1400,700]),'color','w')
    %set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
    
    subplot(3,1,1);
    hold on
    box on
    grid on
    b1 = bar([gc_plot_I(:,iNRMSE), dnn_plot_I(:,iNRMSE)]);
    set(gca,'fontweight','bold','FontSize', 13)
%     plot(gc_plot_I(:,iNRMSE),'-b');
%     plot(dnn_plot_I(:,iNRMSE),'-r');

    % X-Axis
    xlim([0 length(combined_name)+1]);
    xticks([1:length(combined_name)]);
    xticklabels([]);    
    
    % Y-Axis
    ylim([0 ceil(10*nrmse_max)/10]);
    ylabel('NRMSE','fontweight','bold','FontSize', 13)
    title(strcat('Graph Cut vs. DFD-Net Performance Comparision',32, '- Illumination',32,num2str(illum)),'fontweight','bold','FontSize', 16)

    b1(1).FaceColor = [0 0 1];
    b1(2).FaceColor = [1 0 0];
    ax = gca;
    ax.Position = [0.05 0.73 0.92 0.22];    
    
    nmae_max = max(max(gc_data(:,iNMAE)),max(dnn_data(:,iNMAE)));
    subplot(3,1,2);
    hold on
    box on
    grid on
    b2 = bar([gc_plot_I(:,iNMAE), dnn_plot_I(:,iNMAE)]);
    set(gca,'fontweight','bold','FontSize', 13)
%     plot(gc_plot_I(:,iNMAE),'-b');
%     plot(dnn_plot_I(:,iNMAE),'-r');

    % X-Axis
    xlim([0 length(combined_name)+1]);
    xticks([1:length(combined_name)]);
    xticklabels([]);
    
    % Y-Axis
    ylim([0.0 ceil(nmae_max*10)/10])
    ylabel('NMAE','fontweight','bold','FontSize', 13)

    b2(1).FaceColor = [0 0 1];
    b2(2).FaceColor = [1 0 0];
    ax = gca;
    ax.Position = [0.05 0.485 0.92 0.22];  
    
    subplot(3,1,3);
    hold on
    box on
    grid on    
    b3 = bar([gc_plot_I(:,iSSIM), dnn_plot_I(:,iSSIM)]);
    set(gca,'fontweight','bold','FontSize', 13)
%     plot(gc_plot_I(:,iSSIM),'-b');
%     plot(dnn_plot_I(:,iSSIM),'-r')

    % X-Axis
    xlim([0 length(combined_name)+1]);
    xticks([1:length(combined_name)]);
    xticklabels(combined_name);
    xtickangle(90);
    
    % Y-Axis
    ylim([0.0 1.0])
    ylabel('SSIM','fontweight','bold','FontSize', 13)

    b3(1).FaceColor = [0 0 1];
    b3(2).FaceColor = [1 0 0];
    legend('Graph Cut', 'DFD-Net', 'location', 'southoutside', 'orientation', 'horizontal')
    ax = gca;
    ax.Position = [0.05 0.235 0.92 0.22];

    print(plot_num, '-dpng', fullfile(save_dir,strcat('DNN_GC_Comp3_I',num2str(illum),'.png')));

    plot_num = plot_num + 1;

    
end

%% write a table of the comparison between the two

comp_gc_dnn = {};
gc_plot_I2 = {};
dnn_plot_I2 = {};

for edx=1:3
    % get the illumination specific data
    gc_plot_E = gc_data(gc_data(:,2) == edx,:);
    %gc_E_name = gc_data_name(gc_data(:,2) == exp,:);
        
    dnn_plot_E = dnn_data(dnn_data(:,2) == edx,:);
    %dnn_E_name = dnn_data_name(dnn_data(:,2) == ilexplum,:);
    
    for illum=1:3
        gc_met = gc_plot_E(gc_plot_E(:,1) == illum,:);
        dnn_met = dnn_plot_E(dnn_plot_E(:,1) == illum,:);
%         gc_plot_I2{exp,illum} = gc_met;
%         dnn_plot_I2{exp,illum} = dnn_met;
        
        gc_nrmse_mn(edx,illum) = mean(gc_met(:,iNRMSE));
        gc_nmae_mn(edx,illum) = mean(gc_met(:,iNMAE));
        gc_ssim_mn(edx,illum) = mean(gc_met(:,iSSIM));
        
        dnn_nrmse_mn(edx,illum) = mean(dnn_met(:,iNRMSE));
        dnn_nmae_mn(edx,illum) = mean(dnn_met(:,iNMAE));
        dnn_ssim_mn(edx,illum) = mean(dnn_met(:,iSSIM));        
    end
    
    
    gc_ov(edx,1) = mean(gc_plot_E(:,iNRMSE));
    gc_ov(edx,2) = mean(gc_plot_E(:,iNMAE));
    gc_ov(edx,3) = mean(gc_plot_E(:,iSSIM));
    
    dnn_ov(edx,1) = mean(dnn_plot_E(:,iNRMSE));
    dnn_ov(edx,2) = mean(dnn_plot_E(:,iNMAE));
    dnn_ov(edx,3) = mean(dnn_plot_E(:,iSSIM));  
    
end

%% make a table

fprintf('\n\n');
    
caption = strcat('DFD-Net \& Graph Cuts Average Performance Comparison');
write_latex_table_head(caption, 'tbl:dnn_gc_comp1', 'l|l|c|c|c|c|c|c|')
fprintf('\\cline{3-8}\n');
fprintf('\\multicolumn{2}{c|}{\\textbf{}} & \\multicolumn{3}{c|}{\\textbf{Graph Cuts}} & \\multicolumn{3}{c|}{\\textbf{DFD-Net}} \\\\ \\cline{3-8} \n');
fprintf('\\multicolumn{2}{c|}{\\textbf{}} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} & \\textbf{NRMSE} & \\textbf{NMAE} & \\textbf{SSIM} \\\\ \\cline{1-8} \n');

% print out the values
for edx=1:3
    for illum=1:3
        fprintf('\\multicolumn{1}{|l}{Exp %d} &  \\multicolumn{1}{|l|}{Illum %d} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\cline{2-8} \n',edx-1,illum,gc_nrmse_mn(edx,illum),gc_nmae_mn(edx,illum),gc_ssim_mn(edx,illum), dnn_nrmse_mn(edx,illum), dnn_nmae_mn(edx,illum), dnn_ssim_mn(edx,illum));   
    end 
    fprintf('\\multicolumn{1}{|l}{} &  \\multicolumn{1}{|l|}{Overall} & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f & %2.4f \\\\ \\thickhline{2pt} \n',...
        gc_ov(edx,1), gc_ov(edx,2), gc_ov(edx,3), dnn_ov(edx,1), dnn_ov(edx,2), dnn_ov(edx,3));
end
write_latex_table_end;

fprintf('\n\n');

%% compute some other stats

% calculate the means for the graph cuts method
gc_mean = mean(gc_data(:,iNMAE:irun),1);
gc_mean(1,4) = gc_mean(1,4) * 60; %convert to seconds

% calculate the means for the dnn method
dnn_mean = mean(dnn_data(:,iNMAE:irun),1);


fprintf('Graph Cuts Mean: %2.4f, %2.4f, %2.4f, %2.4f\n\n',gc_mean);
fprintf('DfD Net Mean: %2.4f, %2.4f, %2.4f, %2.4f\n\n',dnn_mean);

% calculate the percent differnce between the two

gc_dnn_diff = (gc_mean - dnn_mean)./gc_mean;
fprintf('GC/DNN Difference: %2.2f%%, %2.2f%%, -%2.2f%%, %2.2f%%\n\n',gc_dnn_diff*100);


