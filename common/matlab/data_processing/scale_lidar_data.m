format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% globals

f1 = [5 1];
f2 = [1 5];
f3 = [3 3];
f4 = [5 5];
f5 = [7 7];

% this one gets us to 900x720
x_scale = [2 2 4 4 6 6];
y_scale = [3 6 9 12 15 18];

% this one gets us to 300x240
%x_scale = [1 2 2];
%y_scale = [2 4 6];

% this one does not do any scaling
%x_scale = [];
%y_scale = [];

%% get the directory to scan

startpath = 'D:\IUPUI\Test_Data\';
%img_path = uigetdir(startpath,'Select Folder Containing Scenarios');
img_path = uigetfolder(startpath, 'Select Folders Containing Scenarios');

if(isempty(img_path))
    return;
end

commandwindow;
%listing = dir(img_path);
%listing = listing(3:end);       % remove .  and .. from the listing

%% scan through the directories and get the images
ld_ext = 'bin';
ld_search = strcat('/lidar/*.',ld_ext);
%ld_search = strcat('/lidar/*8bit.',ld_ext);

index = 1;

for idx=1:numel(img_path)
    
    %if(listing(idx).isdir)
        
    %ld_listing = dir(strcat(listing(idx).folder,'/',listing(idx).name,ld_search));
    ld_listing = dir(strcat(img_path{idx},'/',ld_search));

    for jdx=1:numel(ld_listing)
        img_name{index,1} = fullfile(ld_listing(jdx).folder,ld_listing(jdx).name);
        if(strcmp(ld_ext,'bin'))
            img{index,1} = read_binary_lidar_data(fullfile(ld_listing(jdx).folder,ld_listing(jdx).name));
        else
            img{index,1} = imread(fullfile(ld_listing(jdx).folder,ld_listing(jdx).name));           
        end

        % filter the input lidar data
        tmp_img{index,1} = img{index,1};
        img{index,1} = medfilt2(img{index,1}, f1, 'symmetric');
        img{index,1} = medfilt2(img{index,1}, f2, 'symmetric');

        index = index + 1;
    end

    %end 
    
end

%% find the image stats

% get the min
for idx=1:numel(img)  
    tmp_min(idx) = min(img{idx,1}(:));
    tmp_max(idx) = max(img{idx,1}(:));
end

img_min = min(tmp_min);
img_max = max(tmp_max);

fprintf('Min: %2.4f\nMax: %2.4f\n',img_min, img_max);

%%
% figure
% hold on
% plot([img_min:1:img_max],255*(([img_min:1:img_max])-img_min)/(img_max-img_min),'b');
% plot([img_min:1:img_max],(255/img_max)*[img_min:1:img_max],'g');
% plot([img_min:1:img_max],(1/10)*[img_min:1:img_max],'r');


%% scale the images
% img_min = 140;
% img_max = 255;

scale = true;
scale_option = 6;

% use with scale_option = 1,2
% ld_ext = '_8bit_e2.png';          % 8 bit; enhanced contrast
% ld_ext = '_8bit_x3e2.png';        % 8 bit; 2x scale horizontal/vertical; enhanced contrast
% ld_ext = '_8bit_x6e4.png';        % 8 bit; 6x scale in the horizontal, 18x scale in the vertical; enhanced contrast

% use with scale_option = 3
% ld_ext = '_8bit.png';             % 8 bit
% ld_ext = '_8bit_x2.png';          % 8 bit; 2x scale in the horizontal, 3x scale in the vertical; nearest
% ld_ext = '_8bit_x6.png';          % 8 bit; 6x scale in the horizontal, 18x scale in the vertical; nearest

% use with scale_option = 4
% ld_ext = '_8bit_x6d2.png';          % 8 bit; 6x scale in the horizontal, 18x scale in the vertical; nearest decreased res

% use with scale_option = 5
% ld_ext = '_8bit_x6d5.png';          % 8 bit; 6x scale in the horizontal, 18x scale in the vertical; nearest decreased res

% use with scale_option = 6
ld_ext = '_8bit_x6d6.png';          % 8 bit; 6x scale in the horizontal, 18x scale in the vertical; nearest decreased res


figure
for idx=1:numel(img)

    ld_m = img{idx,1};
    
%     % filter the contrast scaled image
%     ld_m = medfilt2(ld_m, f3, 'symmetric');
%     ld_m = medfilt2(ld_m, f4, 'symmetric');
    
    switch scale_option           
        case 1
            % change the scaling based on the max value
            ld_m = (255/img_max)*ld_m;            
        case 2
            % change the contrast scaling of the lidar data
            ld_m = 255*((ld_m)-img_min)/(img_max-img_min);  
        case 3
            ld_m = ld_m/10;
            
        case 4
            ld_m = ld_m/25;
            
        case 5
            ld_m = floor(ld_m/10);
            ld_m = ld_m*10;
            
        case 6
            ld_m = floor(ld_m/5) - 250;
    end
    
    ld_sz = size(ld_m);
    
    % scale up the lidar data 
    if(scale == true)
%         ld_m_up = imresize(ld_m, 2, 'nearest', 'Antialiasing', false);  
%         ld_m_up = medfilt2(ld_m_up, [13 1], 'symmetric');
%         ld_m_up = medfilt2(ld_m_up, [1 13], 'symmetric');          
%         ld_m_s{1} = ld_m_up;

        ld_m_s{1} = ld_m;

        % scale up the lidar data in increments 
        for jdx=1:numel(x_scale)
            % ld_m_s{jdx+1} = imresize(ld_m_s{jdx}, [ld_sz(1)*y_scale(jdx), ld_sz(2)*x_scale(jdx)], 'lanczos3', 'Antialiasing', true);
            ld_m_s{jdx+1} = imresize(ld_m_s{jdx}, [ld_sz(1)*y_scale(jdx), ld_sz(2)*x_scale(jdx)], 'nearest', 'Antialiasing', false);

            % perform final filtering
            %if(jdx==1)
                ld_m_s{jdx+1} = medfilt2(ld_m_s{jdx+1}, [13 1], 'symmetric');                
                ld_m_s{jdx+1} = medfilt2(ld_m_s{jdx+1}, [1 13], 'symmetric');
            %else
                %ld_m_s{jdx+1} = medfilt2(ld_m_s{jdx+1}, f3, 'symmetric');
            %end
        end
        
        ld_m = ld_m_s{end};
    end
    
    image(uint8(ld_m))
    colormap(gray(256));
    axis off;
    drawnow;
    pause(0.5);
    [img_path, img_nm]=fileparts(img_name{idx,1});
    fprintf('Saving: %s\n', strcat(img_path, '/', img_nm, ld_ext));
    
    switch scale_option
        case {1,2,3,4,6}
            imwrite(uint8(ld_m), strcat(img_path, '/', img_nm, ld_ext)); 
        case {5}
            imwrite(uint16(ld_m), strcat(img_path, '/', img_nm, ld_ext));   
    end
end

fprintf('Complete!\n');




