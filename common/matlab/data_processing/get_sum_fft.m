function [img_vals] = get_sum_fft(image_listing)

    for idx=1:numel(image_listing)
        
        img = rgb2gray(imread(fullfile(image_listing(idx).folder, image_listing(idx).name)));
        
        Y = abs(fftshift(fft2(img)))/(size(img,1)*size(img,2));
        img_vals(idx) = sum(Y(:));
    end
    
end