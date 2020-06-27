%function RGB = depth_overlay(img, depth_map)
function depth_overlay(img, depth_map)

    rows = size(depth_map,1);
    cols = size(depth_map,2);
    
    figure();
    set(gcf,'position',([100,100,800,600]),'color','w')
    surf(cols:-1:1,1:rows,depth_map,'CData',img,'FaceColor','texturemap', 'EdgeColor','none');
    view(-160,60);
    hold on
    axis off

    ax = gca;
    ax.Position = [0.03 0.03 0.90 0.92];

end