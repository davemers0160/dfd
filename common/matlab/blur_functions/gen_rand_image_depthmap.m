function [img, dm] = gen_rand_image_depthmap(img_h, img_w, dm_values, color, limits)

    %{circle, polygon, rect}
    circle = limits{1};
    polygon = limits{2};
    rect = limits{3};
    
    % get the random background color
    bg_color = randi([1,numel(color)],1);
    img = cat(3, color{bg_color}(1).*ones(img_h, img_w), color{bg_color}(2).*ones(img_h, img_w), color{bg_color}(3).*ones(img_h, img_w));
    dm = zeros(img_h, img_w, 3);
    
    D = randi([2, numel(dm_values)],1,50);
    D = sort(unique(D));

    for idx=1:numel(D)

        % get the number of shapes for a given depth map value
        N = randi([3,floor(exp((4*(numel(D)-idx))/numel(D))) + 6], 1);

        dm_val = (dm_values(D(idx))/255)*ones(1,3); %, dm_values((idx, idx]/255.0;
        
        for jdx=1:N
            % get the shape type
            T = randi([1,3], 1);
            
            % get the random color for the shape
            C = randi([1,numel(color)],1);
            while(C == bg_color)
                C = randi([1,numel(color)],1);
            end
            
            switch(T)
                case 1
                    X = randi(circle(1,:), 1);
                    Y = randi(circle(2,:), 1);
                    R = randi(circle(3,:), 1);

                    img = insertShape(img, 'FilledCircle', [X, Y, R], 'Color', color{C}, 'Opacity',1, 'SmoothEdges', false);
                    dm = insertShape(dm, 'FilledCircle', [X, Y, R], 'Color', dm_val, 'Opacity',1, 'SmoothEdges', false);

                case 2
                    X = randi(polygon(1,:), 1);
                    Y = randi(polygon(2,:), 1);            
                    P = randi(polygon(3,:), [1,6]);
                    P(1:2:end) = P(1:2:end) + X;
                    P(2:2:end) = P(2:2:end) + Y;
                    P = cat(2, X, Y, P);

                    img = insertShape(img, 'FilledPolygon', P, 'Color', color{C}, 'Opacity',1, 'SmoothEdges', false);
                    dm = insertShape(dm, 'FilledPolygon', P, 'Color', dm_val, 'Opacity',1, 'SmoothEdges', false);
                    
                case 3
                    X = randi(rect(1,:), 1);
                    Y = randi(rect(2,:), 1);                      
                    W = randi(rect(3,:), 1);
                    H = randi(rect(3,:), 1);
                    
                    img = insertShape(img, 'FilledRectangle', [X, Y, W, H], 'Color', color{C}, 'Opacity',1, 'SmoothEdges', false);
                    dm = insertShape(dm, 'FilledRectangle', [X, Y, W, H], 'Color', dm_val, 'Opacity',1, 'SmoothEdges', false);
            end
        end
    end

end