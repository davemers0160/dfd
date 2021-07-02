format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% shape parameters

commandwindow;

% the number of images to generate - not including the intensity variants
num_images = 400;

% depth map values - arrange from lowest to highest with 0 being the lowest
% depthmap_range = [0:1:49];
fg_dm_value = 0;
bg_dm_value = 22;
depthmap_range = [fg_dm_value+1:1:bg_dm_value-1];

min_depthmap = fg_dm_value;
max_depthmap = bg_dm_value;
%num_dm_values = numel(depthmap_range);

%% create all of the image generation parameters

% max number of depth map values for a single image
DM_N = 10;

% initial number of objects at the first depthmap value
num_objects = 50;          

% block dimensions
blk_h = 60;
blk_w = 60;
max_blk_dim = max(blk_h, blk_w);

% build inmage dimensions
img_w = 512 + ceil(blk_w/2);
img_h = 512 + ceil(blk_h/2);
max_dim = max(img_w, img_h);

% make this cropping a mod 16 number
%img_w_range = 20:379;
%img_h_range = 20:379;
img_w_range = ceil(blk_w/2):img_w-1;
img_h_range = ceil(blk_h/2):img_h-1;

% intensity values to simulate different light conditions
%int_values = [0.2, 0.4, 0.6, 0.8, 1.0];
int_values = [1.0];

% x_min,x_max; y_min,y_max; min_r,max_r
rect = [ceil(max_blk_dim/2), ceil(max_blk_dim/0.9)];
circle = [ceil(max_blk_dim/3), ceil(max_blk_dim/1.5)];
polygon = [-ceil(max_blk_dim/2), ceil(max_blk_dim)];
shape_lims = {circle, polygon, rect};

% rect_l = [ceil(max_dim/18), ceil(max_dim/14)];
% circle_l= [ceil(max_dim/18), ceil(max_dim/14)];
% polygon_l = [-ceil(max_dim/11), ceil(max_dim/11)];
% shape_lims_l = {circle_l, polygon_l, rect_l};

% set the probablility of selecting the background depthmap value
prob_bg = 0.325;

% set the probability of selecting the foreground depthmap value
prob_fg = 0.35;

%% start to create the images

% dm_hist = zeros(1, 24);
dm_hist = zeros((max_depthmap - min_depthmap + 1), 1);
dm_hist2 = zeros((max_depthmap - min_depthmap + 1), 1);

tic;
for kdx=0:(num_images-1)
    
    
    % randomly generate DM_N depth values 
    num_dm = DM_N;
        
    % get the probablility that the background depthmap value will be used
    bg_x = rand(1);
    
    % get the probability that the foreground depthmap value will beused 
    fg_x = rand(1);
    
    if(bg_x <= prob_bg)
        num_dm = num_dm - 1;
    end
    
    if(fg_x <= prob_fg)
       num_dm = num_dm - 1;
    end 
    
    if(depthmap_range(end) >= depthmap_range(1))
        D = randi([depthmap_range(1), depthmap_range(end)], 1, num_dm);
    else
        D = randi([depthmap_range(end), depthmap_range(1)], 1, num_dm);
    end    
    
    if(bg_x <= prob_bg)
        D(end+1) = bg_dm_value;
    end
    
    if(fg_x <= prob_fg)
        D(end+1) = fg_dm_value;
    end    
    
    D = sort(unique(D), 'descend');    
    
%     disp(D);
    
    % generate the first depthmap value - use the largest (farthest) value
    dm = (D(1))*ones(img_h, img_w, 1);
   
    for idx=2:numel(D)

        % get the number of shapes for a given depth map value
%         min_N = ceil(exp((3.7*(num_objects-idx))/num_objects));
%         max_N = ceil(exp((4.0*(num_objects-idx))/num_objects));
%         min_N = ceil(num_objects*(exp((4.0*((D(idx)-max_depthmap)))/max_depthmap)));
%         max_N = ceil(num_objects*(exp((1.0*((D(idx)-max_depthmap)))/max_depthmap)));
%         min_N = ceil( ((num_objects)/(1+exp(-0.2*D(idx)+(0.1*num_objects))) ) + 3);
        min_N = ceil( ((max_depthmap)/(1+exp(-0.35*D(idx)+(0.175*max_depthmap))) ) + 3);
%         min_N = floor(D(idx) + 3);

        max_N = ceil(3.0*min_N);
        
        % create the number of objects at a given depth value
        N = randi([min_N, max_N], 1);
        
%         fprintf('%d ', N);
        
        % create the depthmap mask
        dm_mask = zeros(img_h, img_w, 3);
                
%         scale = 1/(max_depthmap + 1 - D(idx));
        
        scale = 1;
        
        for jdx=1:N
            % get the shape type
            T = randi([1,3], 1);

            % get the random color for the shape
            C = [1, 1, 1];
        
            switch(T)
                % filled circle
                case 1
                    X = randi([1,img_w], 1);
                    Y = randi([1,img_h], 1);
                    R = randi(scale * circle(1,:), 1);

                    dm_mask = insertShape(dm_mask, 'FilledCircle', [X, Y, R], 'Color', C, 'Opacity',1, 'SmoothEdges', false);

                % filled polygon
                case 2
                    X = randi([1,img_w], 1);
                    Y = randi([1,img_h], 1);
                    W = randi(scale * rect(1,:), 1);
                    H = randi(scale * rect(1,:), 1);
        
                    S = randi([3, 8], 1);
                    A = 2*pi/S;
                    P = zeros(2, S);
                    
                    for mdx = 0:S-1
                        % a + (b-a).*rand(N,1).
                        angle = mdx * A + ((mdx + 1) * A - mdx * A)*rand(1);

                        if (W / abs(cos(angle)) <= H / abs(sin(angle)))
                            max_radius = ceil(abs(W / cos(angle)));
                        else
                            max_radius = ceil(abs(H / sin(angle)));
                        end

                        radius = randi([floor(max_radius/2), max_radius], 1);
%                         pts.push_back(cv::Point((int32_t)(radius * std::cos((CV_PI / 180.0) * angle)), (int32_t)(radius * std::sin((CV_PI / 180.0) * angle))));
                        P(:, mdx+1) = [floor(radius * cos(angle)); floor(radius * sin(angle))];
                        
                    end
                    
                    P = P(:)';
                    %P = randi(scale * polygon(1,:), [1,6]);
                    P(1:2:end) = P(1:2:end) + X;
                    P(2:2:end) = P(2:2:end) + Y;
                    P = cat(2, X, Y, P);

                    dm_mask = insertShape(dm_mask, 'FilledPolygon', P, 'Color', C, 'Opacity',1, 'SmoothEdges', false);

                % filled rectangle
                case 3
                    A = pi/2 * rand(1);
                    rm = [cos(A), sin(A); -sin(A), cos(A)];
                    
                    X = randi([1,img_w], 1);
                    Y = randi([1,img_h], 1);
                    W = randi(scale * rect(1,:), 1);
                    H = randi(scale * rect(1,:), 1);
                    
                    R = [X, X+W, X+W, X; Y, Y, Y+H, Y+H];
                    R = floor(rm*R);
                    
                    dm_mask = insertShape(dm_mask, 'FilledPolygon', R(:)', 'Color', C, 'Opacity',1, 'SmoothEdges', false);
            end
        end

        dm_mask = dm_mask(:, :, 1);
        
        dm(dm_mask == 1) = D(idx);
        
        bp = 1;
        
    end
    
%     fprintf('\n');

    dm = dm(img_h_range, img_w_range);
    dm_tmp = dm(:);
    
    for idx=1:(max_depthmap - min_depthmap + 1)
        dm_hist(idx,1) = dm_hist(idx,1) + sum(dm_tmp==(idx-1));
    end
       
    for idx=1:(max_depthmap - min_depthmap + 1)
        dm_hist2(idx,1) = dm_hist2(idx,1) + sum(D==(idx-1));
    end
    
    bp = 1;
    
end
toc;

%% histogram plot

min_bin = min_depthmap;
max_bin = max_depthmap;

plot_step = 1;
hist_bin_step = 1;
hist_bins = min_depthmap:hist_bin_step:max_depthmap;

x_lim = [min_depthmap-1,max_depthmap+1];
x = [min_depthmap:hist_bin_step:(max_depthmap)];

y_m = ceil(log10(max(dm_hist(:))));
%y_max = 10^y_m;
y_max = 1.02*max(dm_hist(:))+10;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on
b = bar(x, dm_hist);
%set(gca,'fontweight','bold','FontSize',13, 'yscale', 'log');
set(gca,'fontweight','bold','FontSize',13);

% X-Axis
xlim(x_lim);
%xticks([(x_lim(1)+2):plot_step:(x_lim(2)-2)]);
xticks([min_depthmap:plot_step:max_depthmap]);
xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ylabel('Depth Map Value Count','fontweight','bold')
ylim([1 y_max]);
%ytickformat('%1.2f');

b(1).FaceColor = 'b';
%b(2).FaceColor = 'r';

%title('Depth Map Distribution Comparison', 'fontweight','bold','FontSize',16);
lgd = legend('Ground Truth', 'location','southoutside', 'orientation', 'horizontal');
ax = gca;
ax.Position = [0.06 0.17 0.92 0.78];

plot_num = plot_num + 1;

