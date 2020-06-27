function [img_vals] = get_img_var(image_listing)

    for idx=1:numel(image_listing)
        
        img = double(rgb2gray(imread(fullfile(image_listing(idx).folder, image_listing(idx).name))));
        
        %Y = abs(fftshift(fft2(img)))/(size(img,1)*size(img,2));
        v = var(img(:));
        img_vals(idx) = v;
    end
    
end