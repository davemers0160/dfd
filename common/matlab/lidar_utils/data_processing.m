%%

format long g
format compact
clc
close all
clearvars

full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%v_angle = [...




beam_altitude_angles = [ ...
    16.965853271620, 16.349350684058, 15.776392888927, 15.261876788900, 14.745641815487, 14.147473877370, 13.593996647274, 13.053697446466, ...
    12.527149232740, 11.950180733044, 11.402433080899, 10.885625149691, 10.367098345097, 9.8141940727959, 9.2756137453729, 8.7513573628282, ...
    8.2591866168108, 7.6833640327043, 7.1717127216525, 6.6090681668340, 6.1283565767193, 5.5639931485154, 5.0483311328977, 4.5149074256309, ...
    4.0141423126865, 3.4457681799168, 2.9106255992646, 2.3915258368761, 1.8552373406336, 1.3521803965087, 0.7987031664124, 0.2784574884336, ...
    -0.2584039656040, -0.77864964358279, -1.315511097620, -1.8523725516580, -2.3886610479004, -2.9243765863477, -3.4761349430587, -3.9774730137982, ... 
    -4.5114696788601, -5.02827761006811, -5.577171177803, -6.0922602356260, -6.6394349299760, -7.1516591988229, -7.6799262859336, -8.2397060517764, ...
    -8.7811511681750, -9.30483459292457, -9.843414920348, -10.380276374385, -10.914845997242, -11.448269704509, -11.946742986273, -12.539754304233, ...
    -13.066302517958, -13.6066017187668, -14.14460908839, -14.711264347779, -15.259011999924, -15.773528099952, -16.361955755551, -16.916578901238 ]';

beam_azimuth_angles = [...
        3.25038957, 1.08919277, -1.14247784, -3.27044309, 3.23377380, 1.07257699, -1.07544178, -3.25382732, ...
        3.16673773, 1.03877248, -1.10924629, -3.23721154, 3.18392647, 1.05538826, -1.09205756, -3.17017548, ...
        3.18392647, 1.07257699, -1.07544178, -3.17017548, 3.15012196, 1.03877248, -1.07544178, -3.17017548, ...
        3.16673773, 1.03877248, -1.07544178, -3.17017548, 3.16673773, 1.03877248, -1.02502150, -3.17017548, ...
        3.16673773, 1.07257699, -1.04163727, -3.18679126, 3.20054224, 1.07257699, -1.05882601, -3.15355970, ...
        3.18392647, 1.05538826, -1.05882601, -3.17017548, 3.21715802, 1.08919277, -1.04163727, -3.18679126, ...
        3.20054224, 1.10580854, -1.05882601, -3.22002281, 3.21715802, 1.08919277, -1.05882601, -3.20340703, ...
        3.28419408, 1.08919277, -1.07544178, -3.23721154, 3.30080986, 1.10580854, -1.07544178, -3.28705887]';

h_angle = [-180+(360/4096):360/2048:180-(360/4096)];


lidar_struct = struct('beam_altitude_angles',beam_altitude_angles,'beam_azimuth_angles',beam_azimuth_angles,'h_angle',h_angle);

%     point.x = -r * cos(v_angle) * cos(h_angle);
%     point.y = r * cos(v_angle) * sin(h_angle);
%     point.z = r * sin(v_angle);


%% get images

focus_file = 'D:\IUPUI\Test_Data\Real_World\Garage3\left\exp_30\image_128_30.00_15356234_20180915_121838.png';
defocus_file = 'D:\IUPUI\Test_Data\Real_World\Garage3\left\exp_30\image_133_30.00_15356234_20180915_121838.png';

lidar_file = 'D:\IUPUI\Test_Data\Real_World\Garage3\lidar\lidar_range_00000_20180915_121840.bin';


focus_img = imread(focus_file);
defocus_img = imread(defocus_file);


% crop rectangle as [xmin ymin width height]
left_crop_rect = [914 8 209 54]; %[8 62 914 1123];

gt_img = read_binary_lidar_data(lidar_file);

%gt_img = max(gt_img(:)) - gt_img;

% gt_img2 = medfilt2(gt_img, [3,3]);

% nope
%gt_img_w = wiener2(gt_img,[3,3]);

% gt_img_s = imsharpen(gt_img2,'Radius',2,'Amount',1,'Threshold',0.1);
% 
% 
% gt_img_u = (gt_img - gt_img2);
% gt_img_s2 = gt_img + gt_img_u.*0.9;

% gt_img3 = gt_img2(left_crop(1):left_crop(2),left_crop(3):left_crop(4));

%% plot the two versions of the lidar data

%gt_small = gt_img(left_crop_rect(1):left_crop_rect(2),left_crop_rect(3):left_crop_rect(4));
gt_small = imcrop(gt_img,left_crop_rect);

max_data = max(gt_small(:));

figure(plot_num)
set(gcf,'position',([100,100,600,450]),'color','w')
image(gt_small);
colormap(jet(ceil(max_data)));
title('input');
axis off
ax = gca;
ax.Position = [0.05 0.05 0.90 0.90];
plot_num = plot_num + 1;


% figure(plot_num)
% set(gcf,'position',([100,100,600,450]),'color','w')
% image(gt_img2(left_crop(1):left_crop(2),left_crop(3):left_crop(4)));
% colormap(jet(ceil(max_data)));
% title('median blur');
% axis off
% ax = gca;
% ax.Position = [0.05 0.05 0.90 0.90];
% plot_num = plot_num + 1;

% figure(plot_num)
% set(gcf,'position',([100,100,600,450]),'color','w')
% image(gt_img_u(left_crop(1):left_crop(2),left_crop(3):left_crop(4)));
% colormap(jet(ceil(max(gt_img_u(:)))));
% axis off
% ax = gca;
% ax.Position = [0.05 0.05 0.90 0.90];
% plot_num = plot_num + 1;

% figure(plot_num)
% set(gcf,'position',([100,100,600,450]),'color','w')
% image(gt_img_s2(left_crop(1):left_crop(2),left_crop(3):left_crop(4)));
% colormap(jet(ceil(max_data)));
% title('sharpened');
% axis off
% ax = gca;
% ax.Position = [0.05 0.05 0.90 0.90];
% plot_num = plot_num + 1;

%%
gt_x = zeros(64,2048);
gt_y = zeros(64,2048);
gt_z = zeros(64,2048);

for r=1:64
    for c=1:2048
        gt_x(r,c) = gt_img(r,c)*(cosd(beam_altitude_angles(r))*cosd(h_angle(c)))/25;
        gt_y(r,c) = -gt_img(r,c)*(cosd(beam_altitude_angles(r))*sind(h_angle(c)))/25;
        gt_z(r,c) = gt_img(r,c)*sind(beam_altitude_angles(r))/25;
    end
end

gt_x2 = imcrop(gt_x, left_crop_rect);
%gt_x2 = max(gt_x2(:)) - gt_x2;


gt_xyz(:,:,1) = gt_x;
gt_xyz(:,:,2) = gt_y;
gt_xyz(:,:,3) = gt_z;

RGB = ind2rgb(ceil(gt_x2*40),colormap(jet(ceil(max(gt_x2(:)*40)))));
ptCloud = pointCloud(imcrop(gt_xyz,left_crop_rect),'Color',RGB);

[ptCloudOut,inlierIndices,outlierIndices] = pcdenoise(ptCloud);

%try to get the x plane back out of the ptCloud

gtx3 = [ptCloudOut.Location(:,1);zeros(numel(outlierIndices),1)];
gtx3r = reshape(gtx3,210,55)';

%gt_x2_u = resample(gt_x2,18,1);

% gt_x3_u = resample(gt_x2_u',6,1)';

figure(plot_num)
%scatter3(gt_x(:),gt_y(:),gt_z(:),10,'filled')
pcshow(gt_xyz, 'MarkerSize',10)
title('full');
xlabel('x')
ylabel('y')
plot_num = plot_num + 1;

figure(plot_num)
%scatter3(gt_x(:),gt_y(:),gt_z(:),10,'filled')
pcshow(ptCloudOut, 'MarkerSize',15)
hold on
scatter3(0,0,0,20,'filled');
title('pc')
xlabel('x')
ylabel('y')
plot_num = plot_num + 1;

% figure(plot_num)
% image(gtx3r)
% colormap(jet(ceil(max(gtx3r(:)))));
% title('gtx3r')
% axis off
% plot_num = plot_num + 1;

%%
max_data = max(gt_x2(:));

figure(plot_num)
set(gcf,'position',([100,100,600,450]),'color','w')
image(gt_x2);
colormap(jet(ceil(max_data)));
title('x');
axis off
ax = gca;
ax.Position = [0.05 0.05 0.90 0.90];
plot_num = plot_num + 1;

%%

figure(plot_num)
image(focus_img);
title('focus img');
axis off
plot_num = plot_num + 1;


%% test resize

% f1 = imresize(focus_img,[left_crop(2)-left_crop(1)+1,left_crop(4)-left_crop(3)+1],'bicubic');
g2 = imresize(gt_x2,[size(focus_img,1),size(focus_img,2)],'nearest');

% figure()
% image(f1)
% axis off

figure()
image(g2)
colormap(jet(ceil(max_data)));
title('upsample map');
axis off

%RGB = depthOverlay(f1, gt_img3);

RGB = depthOverlay(focus_img, g2);
title('depth overlay');
